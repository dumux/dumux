// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup CCTpfaDiscretization
 * \brief TODO: Doc me
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_HH

#include <vector>
#include <algorithm>
#include <optional>
#include <utility>
#include <numeric>
#include <cstdint>
#include <memory>
#include <limits>

#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/basegridgeometry.hh>
// #include <dumux/discretization/extrusion.hh>
// #include <dumux/discretization/checkoverlapsize.hh>

#include <dumux/experimental/new_assembly/dumux/common/storage.hh>
#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/scvfseed.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometrycache.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometrylocalview.hh>

namespace Dumux {

// TODO: Extrusion
// TODO: OverlapSize

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Default traits for the `CCTpfaGridGeometry`
 */
template<typename GV>
struct DefaultCCTpfaGridGeometryTraits
: public DefaultMapperTraits<GV>
{
private:
    static constexpr bool isNetworkGrid = int(GV::dimension) < int(GV::dimensionworld);
    static constexpr int numFacesCube = 2*GV::dimension;

public:
    //! Per default, we allow for maximally one hanging node per facet
    static constexpr int maxAdjacentElementLevelDifference = 1;
    //! Per default, we set an upper limit of 8 neighbors on surface/network grids
    static constexpr int maxNumBranchesPerScvf = isNetworkGrid ? 8 : 1;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \todo Doc me
 */
template<typename GV,
         bool cacheGlobally = false,
         Concepts::CCTpfaGridGeometryTraits Traits = DefaultCCTpfaGridGeometryTraits<GV>>
class CCTpfaGridGeometry
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = CCTpfaGridGeometry<GV, cacheGlobally, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;

    using Element = typename GV::template Codim<0>::Entity;
    using GridIndex = typename GV::IndexSet::IndexType;
    using LocalIndex = std::uint_least8_t;

    static constexpr bool isNetworkGrid = int(GV::dimension) < int(GV::dimensionworld);
    static constexpr auto maxNumElementFacets = CCTpfa::Detail::maxNumElementFacets<GV::dimension>;
    static constexpr auto maxNumScvfsPerFacet = CCTpfa::Detail::maxNumScvfsPerFacet<Traits>;
    static constexpr auto maxNumScvfNeighbors = CCTpfa::Detail::maxNumScvfNeighbors<Traits>;

    using FacetScvfIndices = DefaultStorage<GridIndex, maxNumScvfsPerFacet>;
    using ElementFacetsToScvfsMaps = DefaultStorage<FacetScvfIndices, maxNumElementFacets>;
    using Cache = CCTpfa::Detail::GridGeometryCache<ThisType, Traits, cacheGlobally>;

public:
    using ScvfSeed = CCTpfa::ScvfSeed<maxNumScvfNeighbors, GridIndex, LocalIndex>;

    using GridView = GV;
    using LocalView = CCTpfa::GridGeometryLocalView<ThisType, Traits, cacheGlobally>;
    using SubControlVolume = typename LocalView::SubControlVolume;
    using SubControlVolumeFace = typename LocalView::SubControlVolumeFace;
    // using Extrusion = Extrusion_t<Traits>;

    using DiscretizationMethod = DiscretizationMethods::CCTpfa;
    static constexpr DiscretizationMethod discMethod{};
    static constexpr auto maxElementStencilSize = CCTpfa::Detail::maxElementStencilSize<GV::dimension, Traits>;

    explicit CCTpfaGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    {
        // // Check if the overlap size is what we expect
        // if (!CheckOverlapSize<DiscretizationMethod>::isValid(gridView))
        //     DUNE_THROW(Dune::InvalidStateException, "The cctpfa discretization method needs at least an overlap of 1 for parallel computations. "
        //                                              << " Set the parameter \"Grid.Overlap\" in the input file.");

        updateConnectivity_();
        cache_.update(*this);
    }

    friend LocalView localView(const ThisType& gg)
    { return {gg, gg.cache_}; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return facetsToScvfsMap_.size(); }

    //! Return the total number of sub control volumes
    std::size_t numScv() const
    { return facetsToScvfsMap_.size(); }

    //! Return the total number of sub control volume faces
    std::size_t numScvf() const
    {
        return std::accumulate(
            scvfSeeds_.begin(),
            scvfSeeds_.end(),
            std::size_t{0},
            [] (std::size_t current, const ScvfSeed& seed) {
                return current + seed.numNeighbors();
        });
    }

    //! Return the scvf seeds that live on the given facet
    std::ranges::range auto scvfSeeds(const Element& element, unsigned int facetIndex) const
    {
        const auto eIdx = this->elementMapper().index(element);
        return std::views::transform(facetsToScvfsMap_[eIdx][facetIndex], [&] (const auto i) {
            return scvfSeeds_[i];
        });
    }

private:
    void updateConnectivity_()
    {
        resetConnectivity_();
        for (const auto& element : elements(this->gridView()))
            registerElementScvfs_(element, getFacetsToRegister_(element));

        assert(std::ranges::all_of(facetsToScvfsMap_, [] (const auto& facetsToScvfs) {
            return std::ranges::none_of(facetsToScvfs, [] (const auto& scvfIndices) {
                return scvfIndices.empty();
            });
        }));
    }

    void resetConnectivity_()
    {
        facetsToScvfsMap_.resize(this->gridView().size(0), {});
        for (const auto& e : elements(this->gridView()))
        {
            assert(e.subEntities(1) <= maxNumElementFacets);
            facetsToScvfsMap_[this->elementMapper().index(e)].resize(e.subEntities(1));
        }

        scvfSeeds_.clear();
        scvfSeeds_.reserve(this->gridView().size(1));
    }

    auto getFacetsToRegister_(const Element& element) const requires(!isNetworkGrid)
    {
        DefaultStorage<bool, maxNumElementFacets> result(element.subEntities(1), true);
        for (const auto& is : intersections(this->gridView(), element))
            if (is.neighbor())
                result[is.indexInInside()] = element.level() >= is.outside().level();
        return result;
    }

    auto getFacetsToRegister_(const Element& element) const requires(isNetworkGrid)
    {
        DefaultStorage<int, maxNumElementFacets> maxFacetLevel(element.subEntities(1), 0);
        for (const auto& is : intersections(this->gridView(), element))
            if (is.neighbor())
                maxFacetLevel[is.indexInInside()] = std::max(
                    maxFacetLevel[is.indexInInside()],
                    is.outside().level()
                );

        DefaultStorage<bool, maxNumElementFacets> result(element.subEntities(1), false);
        for (const auto [i, level] : enumerate(maxFacetLevel))
            result[i] = element.level() >= level;
        return result;
    }

    template<Concepts::Indexable Storage>
    void registerElementScvfs_(const Element& element,
                               const Storage& facetsToRegister)
    {
        const auto elemIdx = this->elementMapper().index(element);
        for (const auto& is : intersections(this->gridView(), element))
            if (facetsToRegister[is.indexInInside()])
                registerScvf_(is, element, elemIdx);
    }

    void registerScvf_(const typename GV::Intersection& is,
                       const Element& insideElement,
                       const GridIndex insideElemIdx)
    {
        if (!is.neighbor())
            registerBoundaryScvf_(is, insideElemIdx);
        else
            registerInteriorScvf_(is, insideElement, insideElemIdx);
    }

    void registerBoundaryScvf_(const typename GV::Intersection& is,
                               const GridIndex insideElemIdx)
    { pushScvf_(insideElemIdx, is.indexInInside()); }

    void registerInteriorScvf_(const typename GV::Intersection& is,
                               const Element& insideElement,
                               const GridIndex insideElemIdx)
    {
        const auto outsideElement = is.outside();
        const auto outsideElemIdx = this->elementMapper().index(outsideElement);
        auto scvfIndex = getScvfIndex_(insideElemIdx, is.indexInInside());
        if (!scvfIndex)
            addNeighbor_(
                pushScvf_(insideElemIdx, is.indexInInside()),
                outsideElemIdx,
                is.indexInOutside()
            );
        else if (!isAmongNeighbors_(scvfSeeds_[*scvfIndex], outsideElemIdx, is.indexInOutside()))
            addNeighbor_(*scvfIndex, outsideElemIdx, is.indexInOutside());
    }

    std::optional<GridIndex> getScvfIndex_(const GridIndex elemIdx,
                                           const unsigned int facetIdx) const
    {
        const auto it = std::ranges::find_if(
            facetsToScvfsMap_[elemIdx][facetIdx],
            [&] (const ScvfSeed& seed) { return isAmongNeighbors_(seed, elemIdx, facetIdx); },
            [&] (std::integral auto idx) { return scvfSeeds_[idx]; }
        );
        if (it != std::ranges::end(facetsToScvfsMap_[elemIdx][facetIdx]))
            return {*it};
        return std::nullopt;
    }

    std::size_t pushScvf_(const GridIndex insideElemIndex,
                          const LocalIndex insideElemFacetIndex)
    {
        const auto nextFaceIdx = scvfSeeds_.size();
        scvfSeeds_.push_back(ScvfSeed{{insideElemIndex, insideElemFacetIndex}});
        facetsToScvfsMap_[insideElemIndex][insideElemFacetIndex].push_back(nextFaceIdx);
        return nextFaceIdx;
    }

    void addNeighbor_(const GridIndex scvfIndex,
                      const GridIndex elementIndex,
                      const LocalIndex facetIndex)
    {
        assert(scvfIndex < scvfSeeds_.size());
        assert(elementIndex < facetsToScvfsMap_.size());
        assert(facetIndex < facetsToScvfsMap_[elementIndex].size());
        scvfSeeds_[scvfIndex].addOutsideFacet({elementIndex, facetIndex});
        facetsToScvfsMap_[elementIndex][facetIndex].push_back(scvfIndex);
    }

    bool isAmongNeighbors_(const ScvfSeed& scvfSeed,
                           const GridIndex elementIndex,
                           const LocalIndex facetIndex) const
    {
        return std::ranges::any_of(
            scvfSeed.neighboringFacets(),
            [&] (const typename ScvfSeed::Facet& facet) {
                return facet.elementIndex == elementIndex && facet.facetIndex == facetIndex;
        });
    }

    Cache cache_;
    std::vector<ScvfSeed> scvfSeeds_;
    std::vector<ElementFacetsToScvfsMaps> facetsToScvfsMap_;
};

} // end namespace Dumux

#endif
