// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief A gridformat-compatible grid type for CVFE discretizations that
 *        uses DuMux DOF indices directly for VTK point field lookup.
 *
 * Unlike gridformat's LagrangePolynomialGrid this class uses the DuMux
 * dofMapper to identify DOFs, making writing higher-order point fields
 * from a DOF-indexed solution vector efficient (no FE evaluation):
 *
 * \code
 *   IO::GridWriter writer{IO::Format::vtu, gridView, IO::order<2>};
 *   writer.setPointField("velocity", x);   // x indexed by order-2 Lagrange DOF index
 *   writer.write("output");
 * \endcode
 *
 * Parallel correctness: internally the grid builds a sequential VTK point
 * index (0…N_used-1) that covers only DOFs referenced by interior-partition
 * cells, filtering out ghost/overlap DOFs not owned by this process.
 * Each Point object carries both the sequential VTK index (used by gridformat
 * for connectivity) and the DuMux DOF index (used for solution lookup).
 *
 * The VTK Lagrange connectivity inside each element is built using the same
 * dune_to_gfmt_sub_entity ordering as gridformat, so the output is valid
 * VTK Lagrange XML. For methods with multiple DOFs per edge/face (PQ3 and
 * higher), the geometry helper's orientation-consistent dofIndex is used to
 * ensure shared-entity DOFs get the same global index from all elements.
 */
#ifndef DUMUX_IO_CVFE_LAGRANGE_GRID_HH
#define DUMUX_IO_CVFE_LAGRANGE_GRID_HH

#include <config.h>

#if DUMUX_HAVE_GRIDFORMAT

#include <algorithm>
#include <array>
#include <ranges>
#include <tuple>
#include <vector>

#include <gridformat/traits/dune.hpp>
#include <gridformat/grid/traits.hpp>
#include <gridformat/grid/cell_type.hpp>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/discretization/pq2/dofhelper.hh>
#include <dumux/discretization/pq3/dofhelper.hh>

namespace Dumux::IO {

#ifndef DOXYGEN
namespace CVFELagrangeDetail {

// Compute the VTK-ordering sort key for a DOF given its LocalKey.
// Primary sort: highest codim first (corners before edges before faces before interior).
// Secondary: VTK sub-entity slot (from dune_to_gfmt_sub_entity).
// Tertiary: within-slot index, reversed for triangle edge subEntity 1.
//
// Triangle edge subEntity 1 (x=0, Dune direction v0→v2) maps to VTK slot 2 which
// traverses in the opposite direction (v2→v0).  Negating key.index() makes larger
// indices sort first, placing the DOF closer to v2 before the DOF closer to v0.
inline std::tuple<int, int, int> vtkSortKey(const Dune::GeometryType& gt, const auto& key)
{
    using namespace GridFormat::Dune::LagrangeDetail;
    const int codim  = key.codim();
    const int vtkSub = dune_to_gfmt_sub_entity(gt, key.subEntity(), codim);
    const bool rev   = gt.isTriangle() && codim == 1 && key.subEntity() == 1;
    const int within = rev ? -(int(key.index())) : int(key.index());
    return {-codim, vtkSub, within};
}

// Compile-time list mapping polynomial order k to the appropriate DOF helper.
// PQ2LagrangeDofHelper and PQ3LagrangeDofHelper are defined in their respective
// discretization headers and provide the uniform 4-argument dofIndex interface.
template<int k, class GridView>
using PQkDofHelpers = std::conditional_t<k == 2,
                                          Dumux::PQ2LagrangeDofHelper<GridView>,
                                          Dumux::PQ3LagrangeDofHelper<GridView>>;

// Lagrange DOF mapper layout: number of interior Lagrange DOFs per geometry type.
// For a simplex of dimension d: binom(k-1, d); for a cube: (k-1)^d.
template<int k>
int lagrangeLayout(Dune::GeometryType gt, int /*gridDim*/)
{
    const int d = gt.dim();
    if (gt.isCube())
    {
        int r = 1;
        for (int i = 0; i < d; ++i) r *= k - 1;
        return r;
    }
    if (gt.isSimplex())
    {
        if (k - 1 < d) return 0;
        int r = 1;
        for (int i = 0; i < d; ++i) r = r * (k - 1 - i) / (i + 1);
        return r;
    }
    DUNE_THROW(Dune::NotImplemented, "Lagrange layout for geometry type " << gt);
}

} // namespace CVFELagrangeDetail
#endif // DOXYGEN


/*!
 * \ingroup InputOutput
 * \brief A gridformat-compatible Lagrange grid built directly from a GridView
 *        and polynomial order \c k, without requiring a DuMux GridGeometry.
 *
 * VTK point indices equal the Lagrange DOF indices produced by an
 * MCMGMapper with the standard Lagrange layout, so a DOF-indexed solution
 * vector can be written directly as a point field.
 *
 * \tparam GV The Dune grid view type.
 * \tparam k Polynomial order (>= 2).
 */
template<class GV, int k>
class CVFELagrangeGrid
{
    static_assert(k >= 2, "CVFELagrangeGrid requires k >= 2");

    using Element = typename GV::template Codim<0>::Entity;
    using ElementGeometry = typename Element::Geometry;
    using GlobalCoordinate = typename ElementGeometry::GlobalCoordinate;
    using Scalar = typename GV::ctype;
    static constexpr int dim = GV::dimension;

    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
    using FECache = Dune::LagrangeLocalFiniteElementCache<Scalar, Scalar, dim, k>;
    using GeomHelper = CVFELagrangeDetail::PQkDofHelpers<k, GV>;

    struct PointData {
        std::size_t index; //!< Sequential VTK point index (used by gridformat)
        std::size_t dofIndex; //!< Lagrange DOF index (used for field lookup)
        GlobalCoordinate coordinates;
    };

public:
    using GridView = GV;
    using Cell = Element;
    using Point = PointData;
    using Position = GlobalCoordinate;

    explicit CVFELagrangeGrid(const GV& gv)
    : gridView_{gv}
    , dofMapper_{gv, CVFELagrangeDetail::lagrangeLayout<k>}
    { build_(); }

    void update(const GV& gv)
    {
        gridView_ = gv;
        dofMapper_.update(gv);
        build_();
    }

    const GV& gridView() const { return gridView_; }

    std::size_t numberOfPoints() const { return usedDofs_.size(); }
    std::size_t numberOfCells()  const { return numCells_; }

    std::size_t numberOfCellPoints(const Element& e) const
    { return feCache_.get(e.type()).localCoefficients().size(); }

    auto points() const
    {
        return std::views::iota(std::size_t{0}, usedDofs_.size())
            | std::views::transform([this](std::size_t vtkIdx) {
                const auto dofIdx = usedDofs_[vtkIdx];
                return Point{vtkIdx, dofIdx, positions_[dofIdx]};
            });
    }

    auto cells() const
    { return GridFormat::Traits::Cells<GV>::get(gridView_); }

    auto points(const Element& e) const
    {
        const auto eIdx = gridView_.indexSet().index(e);
        const std::size_t n = connectivity_[eIdx].size();
        return std::views::iota(std::size_t{0}, n)
            | std::views::transform([this, eIdx](std::size_t vtk) {
                const auto vtkIdx = connectivity_[eIdx][vtk];
                const auto dofIdx = usedDofs_[vtkIdx];
                return Point{vtkIdx, dofIdx, positions_[dofIdx]};
            });
    }

private:
    void build_()
    {
        const auto& gv = gridView_;
        const std::size_t numAllDofs = dofMapper_.size();

        positions_.assign(numAllDofs, GlobalCoordinate{});
        connectivity_.assign(gv.size(0), {});
        numCells_ = 0;

        std::vector<bool> used(numAllDofs, false);
        const auto& idSet = gv.grid().globalIdSet();

        for (const auto& element : GridFormat::Traits::Cells<GV>::get(gv))
        {
            ++numCells_;
            const auto eIdx = gv.indexSet().index(element);
            const auto geo  = element.geometry();
            const auto gt   = element.type();
            const auto& lc  = feCache_.get(gt).localCoefficients();
            const int n     = static_cast<int>(lc.size());

            std::vector<std::pair<std::tuple<int,int,int>, int>> sortable(n);
            for (int i = 0; i < n; ++i)
                sortable[i] = {CVFELagrangeDetail::vtkSortKey(gt, lc.localKey(i)), i};
            std::sort(sortable.begin(), sortable.end());

            connectivity_[eIdx].resize(n);
            for (int vtk = 0; vtk < n; ++vtk)
            {
                const int i    = sortable[vtk].second;
                const auto& lk = lc.localKey(i);
                const auto g   = GeomHelper::dofIndex(dofMapper_, element, lk, idSet);
                connectivity_[eIdx][vtk] = g;
                positions_[g]            = GeomHelper::dofPosition(geo, lk);
                used[g]                  = true;
            }
        }

        vtk_of_dof_.assign(numAllDofs, std::size_t(-1));
        usedDofs_.clear();
        usedDofs_.reserve(numAllDofs);
        for (std::size_t g = 0; g < numAllDofs; ++g)
            if (used[g]) { vtk_of_dof_[g] = usedDofs_.size(); usedDofs_.push_back(g); }

        for (auto& conn : connectivity_)
            for (auto& idx : conn)
                idx = vtk_of_dof_[idx];
    }

    GV gridView_;
    DofMapper dofMapper_;
    FECache feCache_;
    std::size_t numCells_{0};
    std::vector<GlobalCoordinate>         positions_;
    std::vector<std::vector<std::size_t>> connectivity_;
    std::vector<std::size_t>              usedDofs_;
    std::vector<std::size_t>              vtk_of_dof_;
};


} // namespace Dumux::IO


namespace GridFormat::Traits {

// Traits for CVFELagrangeGrid<GV, k>
template<class GV, int k>
struct Points<Dumux::IO::CVFELagrangeGrid<GV, k>> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>& grid)
    { return grid.points(); }
};

template<class GV, int k>
struct Cells<Dumux::IO::CVFELagrangeGrid<GV, k>> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>& grid)
    { return grid.cells(); }
};

template<class GV, int k>
struct NumberOfPoints<Dumux::IO::CVFELagrangeGrid<GV, k>> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>& grid)
    { return grid.numberOfPoints(); }
};

template<class GV, int k>
struct NumberOfCells<Dumux::IO::CVFELagrangeGrid<GV, k>> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>& grid)
    { return grid.numberOfCells(); }
};

template<class GV, int k>
struct NumberOfCellPoints<Dumux::IO::CVFELagrangeGrid<GV, k>,
                          typename Dumux::IO::CVFELagrangeGrid<GV, k>::Cell> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>& grid,
                    const typename Dumux::IO::CVFELagrangeGrid<GV, k>::Cell& cell)
    { return grid.numberOfCellPoints(cell); }
};

template<class GV, int k>
struct CellPoints<Dumux::IO::CVFELagrangeGrid<GV, k>,
                  typename Dumux::IO::CVFELagrangeGrid<GV, k>::Cell> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>& grid,
                    const typename Dumux::IO::CVFELagrangeGrid<GV, k>::Cell& cell)
    { return grid.points(cell); }
};

template<class GV, int k>
struct CellType<Dumux::IO::CVFELagrangeGrid<GV, k>,
                typename Dumux::IO::CVFELagrangeGrid<GV, k>::Cell> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>&,
                    const typename Dumux::IO::CVFELagrangeGrid<GV, k>::Cell& cell)
    { return GridFormat::Dune::LagrangeDetail::cell_type(cell.type()); }
};

template<class GV, int k>
struct PointCoordinates<Dumux::IO::CVFELagrangeGrid<GV, k>,
                        typename Dumux::IO::CVFELagrangeGrid<GV, k>::Point> {
    static const auto& get(const Dumux::IO::CVFELagrangeGrid<GV, k>&,
                            const typename Dumux::IO::CVFELagrangeGrid<GV, k>::Point& point)
    { return point.coordinates; }
};

template<class GV, int k>
struct PointId<Dumux::IO::CVFELagrangeGrid<GV, k>,
               typename Dumux::IO::CVFELagrangeGrid<GV, k>::Point> {
    static auto get(const Dumux::IO::CVFELagrangeGrid<GV, k>&,
                    const typename Dumux::IO::CVFELagrangeGrid<GV, k>::Point& point)
    { return point.index; }
};

} // namespace GridFormat::Traits


namespace GridFormat::Dune::Traits {

// Expose the underlying GridView so gridformat's Dune::Functions integration
// (set_point_function) can bind local functions to elements.
template<class GV, int k>
struct GridView<Dumux::IO::CVFELagrangeGrid<GV, k>> {
    using type = GV;
    static const auto& get(const Dumux::IO::CVFELagrangeGrid<GV, k>& grid)
    { return grid.gridView(); }
};

} // namespace GridFormat::Dune::Traits

#endif // DUMUX_HAVE_GRIDFORMAT

#endif
