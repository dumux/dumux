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
 * \copydoc Dumux::CCTpfaGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_CACHING_HH
#define DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_CACHING_HH

#include <vector>
#include <cstdint>
#include <cassert>
#include <optional>

#include <dumux/geometry/volume.hh>

#include <dumux/experimental/new_assembly/dumux/common/storage.hh>
#include <dumux/experimental/new_assembly/dumux/geometry/normal.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/fvgridgeometrystorage.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>

#include "detail.hh"

namespace Dumux {

//! Specialization for the case of caching all geometries of the grid view
template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, true, Traits>
: public CCTpfaGridGeometryBase<GV, Traits>
{
    using ThisType = CCTpfaGridGeometry<GV, true, Traits>;
    using ParentType = CCTpfaGridGeometryBase<GV, Traits>;

    static constexpr bool isNetwork = int(GV::dimension) < int(GV::dimensionworld);
    static constexpr auto maxNumFacesPerElement = CCTpfa::Detail::maxNumElementFacets<GV::dimension>
                                                  *CCTpfa::Detail::maxNumFacesPerFacet<Traits>;

    using typename ParentType::FaceSeed;
    using Element = typename GV::template Codim<0>::Entity;
    using ElementGeometry = typename Element::Geometry;
    using Coordinate = typename ElementGeometry::GlobalCoordinate;

    using LocalIndex = std::uint_least8_t;
    using ElementScvfs = Dumux::DefaultStorage<std::size_t, maxNumFacesPerElement>;
    using ScvfNeighbors = Dumux::DefaultStorage<std::size_t, CCTpfa::Detail::maxNumFaceNeighbors<Traits>>;

    using GeometryStorage = FVGridGeometryStorage<
        CCTpfa::Detail::Face<GV>,
        CCTpfa::Detail::ScvWithId<GV, ThisType>,
        CCTpfa::Detail::ScvfWithId<GV, ThisType>
    >;

    struct ScvfSeed
    {
        const FaceSeed& faceSeed;
        LocalIndex idxInNeighbors;
    };

public:
    using SubControlVolume = typename GeometryStorage::SubControlVolume;
    using SubControlVolumeFace = typename GeometryStorage::SubControlVolumeFace;

    class LocalView;

    explicit CCTpfaGridGeometry(const GV& gridView)
    : ParentType(gridView)
    { update_(); }

    friend LocalView localView(const ThisType& gg)
    { return {gg}; }

private:
    friend class LocalView;

    void update_()
    {
        reserve_();
        makeGeometries_();
    }

    void reserve_()
    {
        storage_.reserve({
            .numFaces = this->numFaces(),
            .numScvfs = this->numScvf(),
            .numScvs = this->numScv()
        });
        elementScvfs_.resize(this->gridView().size(0));
        scvfToSeedMap_.reserve(this->numScvf());
    }

    void makeGeometries_()
    {
        for (typename GV::IndexSet::IndexType eIdx = 0; eIdx < this->numScv(); ++eIdx)
        {
            const auto& element = this->element(eIdx);
            const auto& eg = element.geometry();
            storage_.pushScv({eIdx, convexPolytopeVolume(eg), eg.center()});
            for (unsigned int facetIdx = 0; facetIdx < element.subEntities(1); ++facetIdx)
                for (const auto& faceSeed : this->faceSeeds_(eIdx, facetIdx))
                    if (faceSeed.insideFacet().elementIndex == eIdx &&
                        faceSeed.insideFacet().facetIndex == facetIdx)
                        pushScvfs_(faceSeed);
        }
    }

    void pushScvfs_(const FaceSeed& faceSeed) requires(isNetwork)
    {
        LocalIndex nIdx = 0;
        const auto& face = pushFace_(faceSeed);
        const auto& insideFacet = faceSeed.insideFacet();
        pushScvf_(faceSeed, face, getNormal_(insideFacet), nIdx++, insideFacet.elementIndex);
        for (const auto& outsideFacet : faceSeed.outsideFacets())
            pushScvf_(faceSeed, face, getNormal_(outsideFacet), nIdx++, outsideFacet.elementIndex);
    }

    void pushScvfs_(const FaceSeed& faceSeed) requires(!isNetwork)
    {
        LocalIndex nIdx = 0;
        const auto& face = pushFace_(faceSeed);
        const auto n = getNormal_(faceSeed.insideFacet());
        pushScvf_(faceSeed, face, Coordinate{n}, nIdx++, faceSeed.insideFacet().elementIndex);
        for ([[maybe_unused]] const auto& outsideFacet : faceSeed.outsideFacets())
        {
            auto outsideN = n;
            outsideN *= -1.0;
            pushScvf_(faceSeed, face, std::move(outsideN), nIdx++, outsideFacet.elementIndex);
        }
    }

    const typename GeometryStorage::Face& pushFace_(const FaceSeed& faceSeed)
    {
        const auto& element = this->element(faceSeed.insideFacet().elementIndex);
        const auto& faceGeo = element.template subEntity<1>(faceSeed.insideFacet().facetIndex).geometry();
        const auto faceIdx = storage_.pushFace({convexPolytopeVolume(faceGeo), faceGeo.center()});
        return storage_.face(faceIdx);
    }

    void pushScvf_(const FaceSeed& seed,
                   const typename GeometryStorage::Face& face,
                   Coordinate&& normal,
                   LocalIndex neighborIdx,
                   std::size_t insideElementIndex)
    {
        const std::size_t curScvfId = storage_.pushScvf({face, std::move(normal)});
        storage_.scvf(curScvfId).id = curScvfId;
        scvfToSeedMap_.push_back({seed, neighborIdx});
        elementScvfs_[insideElementIndex].push_back(curScvfId);
    }

    Coordinate getNormal_(const typename FaceSeed::Facet& facet) const
    { return getNormal(this->element(facet.elementIndex).geometry(), facet.facetIndex); }

    GeometryStorage storage_;
    std::vector<ElementScvfs> elementScvfs_;
    std::vector<ScvfSeed> scvfToSeedMap_;
};


template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, true, Traits>::LocalView
{
    using IndexType = typename GV::IndexSet::IndexType;

public:
    using GridGeometry = ThisType;
    using Element = typename GridGeometry::Element;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    LocalView(const GridGeometry& gg)
    : gridGeometry_(gg)
    {}

    const GridGeometry& gridGeometry() const
    { return gridGeometry_; }

    void bind(const Element& e) &
    {
        bindElement(e);
        findNeighborScvIndices_();
        findNeighborScvfIndices_();
    }

    void bindElement(const Element& e) &
    {
        clear_();
        element_ = e;
        eIdx_ = gridGeometry_.elementMapper().index(e);
    }

    LocalView bind(const Element& e) && { this->bind(e); return std::move(*this); }
    LocalView bindElement(const Element& e) && { this->bindElement(e); return std::move(*this); }

    friend std::ranges::range auto scvs(const LocalView& lv) { return lv.scvs_(); }
    friend std::ranges::range auto scvfs(const LocalView& lv) { return lv.scvfs_(); }
    friend std::ranges::range auto neighborScvs(const LocalView& lv) { return lv.neighborScvs_(); }
    friend std::ranges::range auto neighborScvfs(const LocalView& lv) { return lv.neighborScvfs_(); }

    const Element& element() const
    { return *element_; }

    auto geometry(const SubControlVolume& scv) const
    { return gridGeometry_.element(scv.dofIndex()).geometry(); }

    auto geometry(const SubControlVolumeFace& scvf) const
    {
        const auto& facet = gridGeometry_.scvfToSeedMap_[scvf.id].faceSeed.insideFacet();
        return gridGeometry_.element(facet.elementIndex).template subEntity<1>(facet.facetIndex).geometry();
    }

    bool onBoundary(const SubControlVolumeFace& scvf) const
    {
        const auto& [seed, idxInNeighbors] = gridGeometry_.scvfToSeedMap_[scvf.id];
        return seed.numOutsideNeighbors() == 0;
    }

    const SubControlVolume& insideScv(const SubControlVolumeFace& scvf) const
    {
        const auto& [seed, idxInNeighbors] = gridGeometry_.scvfToSeedMap_[scvf.id];
        return gridGeometry_.storage_.scv(seed.facet(idxInNeighbors).elementIndex);
    }

    const SubControlVolume& outsideScv(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto& [seed, idxInNeighbors] = gridGeometry_.scvfToSeedMap_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < seedPtr->numNeighbors());
        return gridGeometry_.storage_.scv(seed.facet(i).elementIndex);
    }

    const SubControlVolumeFace& flipScvf(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto& [seed, idxInNeighbors] = gridGeometry_.scvfToSeedMap_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < seedPtr->numNeighbors());
        for (const auto& scvfId : gridGeometry_.elementScvfs_[seed.facet(i).elementIndex])
            if (&gridGeometry_.scvfToSeedMap_[scvfId].faceSeed == &seed)
                return gridGeometry_.storage_.scvf(scvfId);
        DUNE_THROW(Dune::InvalidStateException, "Could not find flip scvf");
    }

private:
    void clear_()
    {
        element_.reset();
        neighborScvIndices_.clear();
        neighborScvfIndices_.clear();
    }

    void findNeighborScvIndices_()
    {
        std::ranges::for_each(gridGeometry_.elementScvfs_[eIdx_], [&] (std::integral auto idx) {
            const auto& faceSeed = gridGeometry_.scvfToSeedMap_[idx].faceSeed;
            std::ranges::for_each(faceSeed.facets(), [&] (const auto& facet) {
                if (facet.elementIndex != eIdx_)
                    neighborScvIndices_.push_back(facet.elementIndex);
            });
        });
    }

    void findNeighborScvfIndices_()
    {
        std::ranges::for_each(neighborScvIndices_, [&] (std::integral auto scvIdx) {
            std::ranges::for_each(gridGeometry_.elementScvfs_[scvIdx], [&] (std::integral auto scvfIdx) {
                const auto& faceSeed = gridGeometry_.scvfToSeedMap_[scvfIdx].faceSeed;
                if (std::ranges::any_of(faceSeed.facets(), [&] (const auto& facet) {
                    return facet.elementIndex == eIdx_;
                }))
                    neighborScvfIndices_.push_back(scvfIdx);
            });
        });
    }

    std::ranges::range auto scvs_() const
    { return std::views::single(gridGeometry_.storage_.scv(eIdx_)); }

    std::ranges::range auto scvfs_() const
    {
        return std::views::transform(
            gridGeometry_.elementScvfs_[eIdx_],
            [&] (std::integral auto scvfIdx) -> const SubControlVolumeFace& {
                return gridGeometry_.storage_.scvf(scvfIdx);
        });
    }

    std::ranges::range auto neighborScvs_() const
    {
        return std::views::transform(
            neighborScvIndices_,
            [&] (std::integral auto scvIdx) -> const SubControlVolume& {
                return gridGeometry_.storage_.scv(scvIdx);
        });
    }

    std::ranges::range auto neighborScvfs_() const
    {
        return std::views::transform(
            neighborScvfIndices_,
            [&] (std::integral auto scvfIdx) -> const SubControlVolumeFace& {
                return gridGeometry_.storage_.scvf(scvfIdx);
        });
    }

    const GridGeometry& gridGeometry_;
    std::optional<Element> element_;
    IndexType eIdx_;
    DefaultStorage<IndexType, GridGeometry::maxElementStencilSize> neighborScvIndices_;
    DefaultStorage<IndexType, GridGeometry::maxElementStencilSize> neighborScvfIndices_;
};

} // end namespace Dumux

#endif
