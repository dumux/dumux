// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Dispatches dynamic (post-refinement) load balancing to the grid-specific API.
 *
 * DuMux itself never calls loadBalance() again after the initial one done at grid setup
 * (dumux/adaptive/adapt.hh only ever calls preAdapt/adapt/postAdapt) -- there is no dynamic
 * rebalancing anywhere in the tree today, for any grid. This is new. Both grid backends used
 * here were checked directly against their vendored sources for this plan and already expose
 * what's needed at the Dune level:
 *
 *  - Dune::ALUGrid: `grid.loadBalance(weights, dataHandle)` -- a single call that both
 *    computes a new partition (steered by a per-macro-element `weights` functor) and
 *    migrates the payload carried by `dataHandle`. Demonstrated end-to-end in
 *    dune-alugrid/examples/loadbalance/{adaptation.hh,loadbalance_simple.hh}; the
 *    `LeafCountWeights` functor below is a direct port of that example's
 *    `SimpleLoadBalanceWeights` (weight = number of leaf descendants under a macro element,
 *    i.e. "how much refined work lives here" -- the right proxy for post-refinement
 *    imbalance).
 *  - Dune::UGGrid: unlike ALUGrid, UGGrid's plain `loadBalance()` (the RCB-based algorithm
 *    behind the single-argument `loadBalance(dataHandle)` overload) is a one-shot algorithm
 *    for distributing an initially-serial grid and unconditionally throws
 *    `Dune::NotImplemented` if called again on a grid that is already parallel (see
 *    `BalanceGridRCB` in dune-uggrid/dune/uggrid/parallel/dddif/lbrcb.cc -- confirmed by
 *    actually hitting that throw when this file first called the wrong overload). The only
 *    UGGrid API that *can* redistribute an already-distributed grid is
 *    `loadBalance(targetProcessors, fromLevel, dataHandle)`, which needs the caller to supply
 *    a target rank for every interior leaf element; `computeSpaceFillingCurveTargetProcessors`
 *    below builds that vector via an unweighted Morton-order element-count balance (gathers
 *    all interior element centroids -- small, test-scale data -- so every rank can compute the
 *    identical global assignment independently, no data-handle round trip needed for the
 *    assignment itself, only for the payload `dataHandle` carries).
 *
 * Both branches are only reachable if the corresponding Dune module was found by CMake
 * (HAVE_DUNE_ALUGRID / HAVE_DUNE_UGGRID), matching how the rest of DuMux gates grid-manager
 * code (see dumux/io/grid/gridmanager_{alu,ug}.hh).
 */
#ifndef DUMUX_BOUSSINESQ_ADAPTIVE_LOADBALANCER_HH
#define DUMUX_BOUSSINESQ_ADAPTIVE_LOADBALANCER_HH

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>
#include <set>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/common/gridenums.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#include <dune/grid/common/mcmgmapper.hh>
#endif

namespace Dumux::BoussinesqAdaptive {

#if HAVE_DUNE_ALUGRID
/*!
 * \brief Per-macro-element weight for ALUGrid's weighted load balancing: the number of leaf
 *        descendants, i.e. how much refined work currently lives under this macro element.
 *        Direct port of dune-alugrid/examples/loadbalance/loadbalance_simple.hh's
 *        SimpleLoadBalanceWeights.
 */
template<class Grid>
struct LeafCountWeights
{
    using Element = typename Grid::template Codim<0>::Entity;

    explicit LeafCountWeights(const Grid& grid) : grid_(grid) {}

    long int operator()(const Element& element) const
    {
        const int maxLevel = grid_.maxLevel();
        long int leafElements = 1;
        for (auto it = element.hbegin(maxLevel); it != element.hend(maxLevel); ++it)
            if (it->isLeaf())
                ++leafElements;
        return leafElements;
    }

private:
    const Grid& grid_;
};
#endif

#if HAVE_DUNE_ALUGRID
/*!
 * \brief User-defined ALUGrid partitioning into vertical strips: rank r owns one contiguous
 *        x-interval, each spanning the full domain height ("slices from top to bottom").
 *
 * Motivation (Rayleigh-Bénard fingering): the refined boundary layer spans the full domain
 * width and the fingers grow vertically. The default space-filling-curve partition produces
 * compact 2D patches, so after refinement the rank(s) owning the top boundary layer carry
 * most of the load and every finger crosses partition boundaries repeatedly as it descends.
 * Vertical strips give every rank an equal share of the boundary layer, and a finger stays
 * on its rank for its whole descent -- partition boundaries are only crossed by the (weak)
 * horizontal coupling.
 *
 * The strip cut positions are weighted: each macro element counts with its number of leaf
 * descendants (same load proxy as LeafCountWeights), so on a refined grid the strips have
 * equal *work*, not equal width. Computed identically and deterministically on every rank
 * (one allgatherv of (x, weight) pairs of all interior macro elements -- macro-grid-sized,
 * small), so no result communication is needed. Usable both for the initial partition
 * (weights all 1 on an unrefined grid => equal-width strips) and for dynamic rebalancing
 * via grid.repartition(destinations, dataHandle), where it preserves the strip layout that
 * ALUGrid's own SFC-based loadBalance() would otherwise destroy on the first rebalance.
 */
template<class Grid>
class VerticalStripDestinations
{
public:
    explicit VerticalStripDestinations(const Grid& grid)
    : numStrips_(grid.comm().size())
    {
        LeafCountWeights<Grid> weights(grid);
        const auto macroView = grid.levelGridView(0);
        const auto& comm = grid.comm();
        const int numRanks = comm.size();

        // flattened (xCenter, leafCountWeight) pairs of this rank's interior macro elements
        std::vector<double> local;
        for (const auto& element : elements(macroView, Dune::Partitions::interior))
        {
            local.push_back(element.geometry().center()[0]);
            local.push_back(static_cast<double>(weights(element)));
        }

        const int localCount = static_cast<int>(local.size());
        std::vector<int> counts(numRanks);
        comm.allgather(&localCount, 1, counts.data());
        std::vector<int> displ(numRanks, 0);
        for (int r = 1; r < numRanks; ++r)
            displ[r] = displ[r-1] + counts[r-1];

        std::vector<double> global(displ.back() + counts.back());
        comm.allgatherv(local.data(), localCount, global.data(), counts.data(), displ.data());

        // deterministic global x-order (tie-break on gather position, identical on all ranks)
        const std::size_t n = global.size()/2;
        std::vector<std::size_t> order(n);
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
        { return global[2*a] != global[2*b] ? global[2*a] < global[2*b] : a < b; });

        double totalWeight = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            totalWeight += global[2*i + 1];

        // walk the x-sorted elements, cutting whenever the accumulated weight crosses the
        // next equal-share threshold; the crossing element itself stays left of its cut
        double acc = 0.0;
        int strip = 1;
        for (std::size_t k = 0; k < n && strip < numStrips_; ++k)
        {
            acc += global[2*order[k] + 1];
            while (strip < numStrips_ && acc >= totalWeight*strip/numStrips_)
            {
                cuts_.push_back(global[2*order[k]]);
                ++strip;
            }
        }
    }

    //! target rank of a (macro) element: the number of cuts strictly below its x-center --
    //! consistent with the cut construction above, where the threshold-crossing element
    //! itself belongs to the strip left of the cut it created
    template<class Element>
    int operator()(const Element& element) const
    {
        const double x = element.geometry().center()[0];
        const auto stripIdx = std::lower_bound(cuts_.begin(), cuts_.end(), x) - cuts_.begin();
        return static_cast<int>(stripIdx);
    }

    //! ALUGrid queries which ranks this rank may receive elements from during
    //! repartitioning; returning false lets ALUGrid determine that itself via a global
    //! communication -- required for correctness on the first transition from the
    //! default SFC layout to strips, where elements can arrive from any rank (same
    //! default as dune-alugrid's own loadbalance_simple.hh example)
    bool importRanks(std::set<int>& /*ranks*/) const { return false; }

private:
    int numStrips_;
    std::vector<double> cuts_;
};
#endif

//! Primary template: only the ALUGrid/UGGrid specializations below are usable.
template<class Grid>
struct GridLoadBalancer;

#if HAVE_DUNE_ALUGRID
template<int dim, int dimworld, Dune::ALUGridElementType elType,
         Dune::ALUGridRefinementType refType, class Comm>
struct GridLoadBalancer<Dune::ALUGrid<dim, dimworld, elType, refType, Comm>>
{
    using Grid = Dune::ALUGrid<dim, dimworld, elType, refType, Comm>;

    //! Rebalance the grid, migrating whatever `dataHandle` carries along with it.
    //! \param verticalStrips use the weighted vertical-strip partitioner above (and keep the
    //!        strip layout across rebalances) instead of ALUGrid's default SFC partitioning
    //! \return true if the partitioning actually changed (matches loadBalance()'s own return)
    template<class DataHandle>
    static bool apply(Grid& grid, DataHandle& dataHandle, bool verticalStrips = false)
    {
        if (grid.comm().size() <= 1)
            return false;

        if (verticalStrips)
        {
            VerticalStripDestinations<Grid> destinations(grid);
            return grid.repartition(destinations, dataHandle);
        }

        LeafCountWeights<Grid> weights(grid);
        return grid.loadBalance(weights, dataHandle);
    }
};
#endif

#if HAVE_DUNE_UGGRID
/*!
 * \brief Compute a target rank for every interior leaf element of an already-distributed
 *        UGGrid, balancing element counts across ranks along a Z-order (Morton) curve.
 *
 * UGGrid's own internal partitioner (`grid.loadBalance()`, used by the plain
 * `loadBalance(dataHandle)` overload) is a one-shot algorithm for distributing an initially
 * serial grid and explicitly refuses to redistribute a grid that is already parallel
 * (`BalanceGridRCB` throws `Dune::NotImplemented` for that case, see
 * dune-uggrid/dune/uggrid/parallel/dddif/lbrcb.cc). The only UGGrid API that *can*
 * redistribute is `loadBalance(targetProcessors, fromLevel, dataHandle)`, which requires the
 * caller to supply the target rank for every interior leaf element -- this function computes
 * that vector.
 *
 * Every rank runs the exact same deterministic computation: gather every rank's interior
 * element centroids (small, test-scale data -- one allgatherv, no per-timestep dependency on
 * a runtime-varying message size class), derive a Morton code per centroid normalized against
 * the combined bounding box of the gathered data (no separate reduction needed, the bounding
 * box falls out of the already-global list), sort by (code, tie-break index), and hand out
 * contiguous equal-sized chunks of that order to ranks 0..P-1. Every rank ends up with the
 * identical global assignment and only needs to read back the slice of it that was its own
 * contribution -- no result needs to be communicated back.
 */
template<class Grid>
std::vector<typename Grid::Rank> computeSpaceFillingCurveTargetProcessors(const Grid& grid)
{
    static constexpr int dim = Grid::dimension;
    using GridView = typename Grid::LeafGridView;
    using Rank = typename Grid::Rank;

    const GridView gridView = grid.leafGridView();
    const auto& comm = gridView.comm();
    const int numRanks = comm.size();
    const int rank = comm.rank();

    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper(gridView, Dune::mcmgElementLayout());

    // 1. this rank's own interior elements: mapper index + centroid
    std::vector<int> localIndices;
    std::vector<double> localCenters;
    for (const auto& element : elements(gridView, Dune::Partitions::interior))
    {
        localIndices.push_back(static_cast<int>(elementMapper.index(element)));
        const auto center = element.geometry().center();
        for (int d = 0; d < dim; ++d)
            localCenters.push_back(center[d]);
    }
    const int localCount = static_cast<int>(localIndices.size());

    // 2. counts/displacements (in element units) across all ranks
    std::vector<int> counts(numRanks);
    comm.allgather(&localCount, 1, counts.data());

    std::vector<int> displ(numRanks, 0);
    for (int r = 1; r < numRanks; ++r)
        displ[r] = displ[r-1] + counts[r-1];
    const int totalCount = displ.back() + counts.back();

    // 3. gather all centroids (flattened, dim components per element)
    std::vector<int> countsScaled(numRanks), displScaled(numRanks);
    for (int r = 0; r < numRanks; ++r)
    {
        countsScaled[r] = counts[r]*dim;
        displScaled[r] = displ[r]*dim;
    }
    std::vector<double> globalCenters(totalCount*dim);
    comm.allgatherv(localCenters.data(), localCount*dim,
                     globalCenters.data(), countsScaled.data(), displScaled.data());

    // 4. Morton code per global element, quantized against the bounding box of the
    //    gathered centroids themselves (identical on every rank, no extra reduction needed)
    std::array<double, dim> lo, hi;
    for (int d = 0; d < dim; ++d)
    {
        lo[d] = globalCenters[d];
        hi[d] = globalCenters[d];
    }
    for (int i = 1; i < totalCount; ++i)
        for (int d = 0; d < dim; ++d)
        {
            using std::min; using std::max;
            lo[d] = min(lo[d], globalCenters[i*dim + d]);
            hi[d] = max(hi[d], globalCenters[i*dim + d]);
        }

    static constexpr int bitsPerDim = 63/dim;
    static constexpr std::uint32_t quantMax = (bitsPerDim < 32) ? ((std::uint32_t{1} << bitsPerDim) - 1)
                                                                 : std::numeric_limits<std::uint32_t>::max();

    std::vector<std::uint64_t> mortonCode(totalCount, 0);
    for (int i = 0; i < totalCount; ++i)
    {
        std::array<std::uint32_t, dim> q;
        for (int d = 0; d < dim; ++d)
        {
            const double extent = hi[d] - lo[d];
            const double normalized = (extent > 0.0) ? (globalCenters[i*dim + d] - lo[d])/extent : 0.0;
            q[d] = static_cast<std::uint32_t>(normalized*quantMax + 0.5);
        }

        std::uint64_t code = 0;
        for (int b = 0; b < bitsPerDim; ++b)
            for (int d = 0; d < dim; ++d)
                if (q[d] & (std::uint32_t{1} << b))
                    code |= std::uint64_t{1} << (b*dim + d);
        mortonCode[i] = code;
    }

    // 5. deterministic global order (tie-break on index so equal codes don't depend on sort stability)
    std::vector<int> order(totalCount);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b)
    { return mortonCode[a] != mortonCode[b] ? mortonCode[a] < mortonCode[b] : a < b; });

    // 6. hand out contiguous, as-equal-as-possible chunks of that order to ranks 0..P-1
    std::vector<Rank> targetOfGlobalIndex(totalCount);
    {
        const int base = totalCount/numRanks;
        const int remainder = totalCount % numRanks;
        int pos = 0;
        for (int r = 0; r < numRanks; ++r)
        {
            const int chunkSize = base + (r < remainder ? 1 : 0);
            for (int k = 0; k < chunkSize; ++k)
                targetOfGlobalIndex[order[pos++]] = static_cast<Rank>(r);
        }
    }

    // 7. every rank reads back only the slice of the global assignment it contributed --
    //    no communication of the result is needed, the computation above was identical
    //    (and deterministic) on every rank.
    std::vector<Rank> targetProcessors(gridView.size(0), static_cast<Rank>(rank));
    for (int k = 0; k < localCount; ++k)
        targetProcessors[localIndices[k]] = targetOfGlobalIndex[displ[rank] + k];

    return targetProcessors;
}

template<int dim>
struct GridLoadBalancer<Dune::UGGrid<dim>>
{
    using Grid = Dune::UGGrid<dim>;

    //! \note verticalStrips is accepted for interface parity with the ALUGrid
    //!       specialization but not implemented for UGGrid (Morton-order fallback)
    template<class DataHandle>
    static bool apply(Grid& grid, DataHandle& dataHandle, bool /*verticalStrips*/ = false)
    {
        if (grid.comm().size() <= 1)
            return false;

        const auto targetProcessors = computeSpaceFillingCurveTargetProcessors(grid);
        return grid.loadBalance(targetProcessors, /*fromLevel=*/0, dataHandle);
    }
};
#endif

} // end namespace Dumux::BoussinesqAdaptive

#endif
