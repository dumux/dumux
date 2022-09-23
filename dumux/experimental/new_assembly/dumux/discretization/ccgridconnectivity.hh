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
 * \ingroup CCDiscretization
 * \copydoc Dumux::CCGridConnectivity
 */
#ifndef DUMUX_DISCRETIZATION_CC_GRID_CONNECTIVITY_HH
#define DUMUX_DISCRETIZATION_CC_GRID_CONNECTIVITY_HH

#include <cstdint>
#include <cassert>
#include <ranges>

#include <dumux/common/enumerate.hh>

#include <dumux/experimental/new_assembly/dumux/common/size.hh>
#include <dumux/experimental/new_assembly/dumux/common/storage.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/ccfaceseed.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<std::size_t dim>
inline constexpr std::size_t maxNumElementFacets = 2*dim; // cubes

} // namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup CCDiscretization
 * \brief Class to store the connectivity of the faces of a grid view
 *        in the context of cell-centered discretization schemes.
 */
template<typename GridView,
         Concepts::Size auto maxNumFaceNeighbors = Size::dynamic,
         Concepts::Size auto maxNumFacesPerFacet = Size::dynamic,
         Concepts::Size auto maxNumElementFacets = Detail::maxNumElementFacets<GridView::dimension>>
class CCGridConnectivity
{
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndex = typename GridView::IndexSet::IndexType;
    using LocalIndex = std::uint_least8_t;

    using FacetFaces = DefaultStorage<GridIndex, maxNumFacesPerFacet>;
    using ElementFacetsToFacesMap = DefaultStorage<FacetFaces, maxNumElementFacets>;

    static constexpr bool isNetworkGrid = int(GridView::dimension) < int(GridView::dimensionworld);

public:
    using FaceSeed = Dumux::CCFaceSeed<maxNumFaceNeighbors, GridIndex, LocalIndex>;
    using Facet = typename FaceSeed::Facet;

    template<Concepts::ElementMapper<GridView> Mapper>
    explicit CCGridConnectivity(const GridView& gridView, const Mapper& mapper)
    : gridView_(gridView)
    { update_(mapper); }

    //! Return the number of faces of the discretization
    std::size_t numFaces() const
    { return faceSeeds_.size(); }

    //! Return all face seeds of the discretization
    const auto& faceSeeds() const
    { return faceSeeds_; }

    //! Return seeds of the faces associated with the given facet
    std::ranges::range auto faceSeeds(const Facet& facet) const
    {
        const auto& faceIndices = elementFacetsToFacesMap_[facet.elementIndex][facet.facetIndex];
        return std::views::transform(faceIndices, [&] (const auto i) -> const FaceSeed& {
            return faceSeeds_[i];
        });
    }

    //! Update the connectivity for the given grid view & mapper
    template<Concepts::ElementMapper<GridView> Mapper>
    void update(const GridView& gridView, const Mapper& mapper)
    {
        gridView_ = gridView;
        update_(mapper);
    }

private:
    template<Concepts::ElementMapper<GridView> Mapper>
    void update_(const Mapper& mapper)
    {
        reset_(mapper);

        // create markers once to avoid repeated allocation
        DefaultStorage<int, maxNumElementFacets> facetLevels;
        DefaultStorage<bool, maxNumElementFacets> facetsToRegister;
        for (const auto& element : elements(gridView_))
        {
            markFacetsToRegister_(element, facetLevels, facetsToRegister);
            registerElementFaces_(element, mapper, facetsToRegister);
        }

        assert(isValid_());
    }

    bool isValid_() const
    {
        return std::ranges::all_of(elementFacetsToFacesMap_, [] (const auto& facetsToFaces) {
            return std::ranges::none_of(facetsToFaces, [] (const auto& faceIndices) {
                return faceIndices.empty();
            });
        });
    }

    template<Concepts::ElementMapper<GridView> Mapper>
    void reset_(const Mapper& mapper)
    {
        faceSeeds_.clear();
        faceSeeds_.reserve(gridView_.size(1));

        elementFacetsToFacesMap_.resize(gridView_.size(0), {});
        for (const auto& e : elements(gridView_))
        {
            assert(Size::isRepresentableBy(maxNumElementFacets, e.subEntities(1)));
            elementFacetsToFacesMap_[mapper.index(e)].resize(e.subEntities(1));
        }
    }

    void markFacetsToRegister_(const Element& element,
                               DefaultStorage<int, maxNumElementFacets>& levels,
                               DefaultStorage<bool, maxNumElementFacets>& markers) const requires(!isNetworkGrid)
    {
        markers.assign(element.subEntities(1), true);
        for (const auto& is : intersections(gridView_, element))
            if (is.neighbor())
                markers[is.indexInInside()] = element.level() >= is.outside().level();
    }

    void markFacetsToRegister_(const Element& element,
                               DefaultStorage<int, maxNumElementFacets>& maxFacetLevels,
                               DefaultStorage<bool, maxNumElementFacets>& markers) const requires(isNetworkGrid)
    {
        maxFacetLevels.assign(element.subEntities(1), 0);
        for (const auto& is : intersections(gridView_, element))
            if (is.neighbor())
                maxFacetLevels[is.indexInInside()] = std::max(
                    maxFacetLevels[is.indexInInside()],
                    is.outside().level()
                );

        markers.resize(element.subEntities(1));
        for (const auto [i, level] : enumerate(maxFacetLevels))
            markers[i] = element.level() >= level;
    }

    template<Concepts::ElementMapper<GridView> Mapper>
    void registerElementFaces_(const Element& element,
                               const Mapper& mapper,
                               const DefaultStorage<bool, maxNumElementFacets>& facetsToRegister)
    {
        const auto elemIdx = mapper.index(element);
        for (const auto& is : intersections(gridView_, element))
            if (facetsToRegister[is.indexInInside()])
                registerFace_(is, element, mapper, elemIdx);
    }

    template<Concepts::ElementMapper<GridView> Mapper>
    void registerFace_(const typename GridView::Intersection& is,
                       const Element& insideElement,
                       const Mapper& mapper,
                       const GridIndex insideElemIdx)
    {
        if (!is.neighbor())
            registerBoundaryFace_(is, insideElemIdx);
        else
            registerInteriorFace_(is, insideElement, mapper, insideElemIdx);
    }

    void registerBoundaryFace_(const typename GridView::Intersection& is,
                               const GridIndex insideElemIdx)
    { pushFace_(insideElemIdx, is.indexInInside(), is.boundary()); }

    template<Concepts::ElementMapper<GridView> Mapper>
    void registerInteriorFace_(const typename GridView::Intersection& is,
                               const Element& insideElement,
                               const Mapper& mapper,
                               const GridIndex insideElemIdx)
    {
        const auto outsideElement = is.outside();
        const auto outsideElemIdx = mapper.index(outsideElement);
        auto faceIndex = getFaceIndex_(insideElemIdx, is.indexInInside());
        if (!faceIndex)
        {
            faceIndex = pushFace_(insideElemIdx, is.indexInInside(), is.boundary());
            addNeighbor_(*faceIndex, outsideElemIdx, is.indexInOutside());
        }
        else if (!isAmongNeighbors_(faceSeeds_[*faceIndex], outsideElemIdx, is.indexInOutside()))
            addNeighbor_(*faceIndex, outsideElemIdx, is.indexInOutside());
    }

    std::optional<GridIndex> getFaceIndex_(const GridIndex elemIdx,
                                           const unsigned int facetIdx) const
    {
        const auto it = std::ranges::find_if(
            elementFacetsToFacesMap_[elemIdx][facetIdx],
            [&] (const FaceSeed& s) { return isAmongNeighbors_(s, elemIdx, facetIdx); },
            [&] (std::integral auto idx) { return faceSeeds_[idx]; }
        );
        if (it != std::ranges::end(elementFacetsToFacesMap_[elemIdx][facetIdx]))
            return {*it};
        return std::nullopt;
    }

    std::size_t pushFace_(const GridIndex insideElemIndex,
                          const LocalIndex insideElemFacetIndex,
                          bool boundary)
    {
        const auto nextFaceIdx = faceSeeds_.size();
        faceSeeds_.emplace_back(Facet{insideElemIndex, insideElemFacetIndex}, boundary);
        elementFacetsToFacesMap_[insideElemIndex][insideElemFacetIndex].push_back(nextFaceIdx);
        return nextFaceIdx;
    }

    void addNeighbor_(const GridIndex faceIndex,
                      const GridIndex elementIndex,
                      const LocalIndex facetIndex)
    {
        assert(faceIndex < faceSeeds_.size());
        assert(elementIndex < elementFacetsToFacesMap_.size());
        assert(facetIndex < elementFacetsToFacesMap_[elementIndex].size());
        faceSeeds_[faceIndex].addOutsideFacet({elementIndex, facetIndex});
        elementFacetsToFacesMap_[elementIndex][facetIndex].push_back(faceIndex);
    }

    bool isAmongNeighbors_(const FaceSeed& seed,
                           const GridIndex elementIndex,
                           const LocalIndex facetIndex) const
    {
        return std::ranges::any_of(seed.facets(), [&] (const Facet& facet) {
            return facet.elementIndex == elementIndex && facet.facetIndex == facetIndex;
        });
    }

    const GridView& gridView_;
    std::vector<FaceSeed> faceSeeds_;
    std::vector<ElementFacetsToFacesMap> elementFacetsToFacesMap_;
};

} // end namespace Dumux

#endif
