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
#ifndef DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_NOCACHING_HH
#define DUMUX_DISCRETIZATION_CCTPFA_GRID_GEOMETRY_NOCACHING_HH

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

//! Specialization for the case of caching geometries only locally
template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, false, Traits>
: public CCTpfaGridGeometryBase<GV, Traits>
{
    using ParentType = CCTpfaGridGeometryBase<GV, Traits>;

public:
    class LocalView;
    using SubControlVolume = typename LocalView::SubControlVolume;
    using SubControlVolumeFace = typename LocalView::SubControlVolumeFace;

    explicit CCTpfaGridGeometry(const GV& gridView)
    : ParentType(gridView)
    {}

    friend LocalView localView(const CCTpfaGridGeometry& gg)
    { return {gg}; }

private:
    friend class LocalView;
};


template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, false, Traits>::LocalView
{
    using GG = CCTpfaGridGeometry<GV, false, Traits>;
    using ThisType = typename GG::LocalView;

    static constexpr int dim = GV::dimension;
    static constexpr auto maxNumFaceNeighbors = CCTpfa::Detail::maxNumFaceNeighbors<Traits>;
    static constexpr auto maxNumStencilScvfs = CCTpfa::Detail::maxNumStencilScvfs<dim, Traits>;
    static constexpr auto maxNumElementFaces = CCTpfa::Detail::maxNumElementFaces<dim, Traits>;
    static constexpr auto maxElementStencilSize = CCTpfa::Detail::maxElementStencilSize<dim, Traits>;

    struct StorageSizes
    {
        static constexpr auto maxNumScv = maxElementStencilSize;
        static constexpr auto maxNumScvf = maxNumStencilScvfs;
        static constexpr auto maxNumFaces = maxNumElementFaces;
    };

    using LocalFaceNeighbors = DefaultStorage<unsigned int, maxNumFaceNeighbors>;

    using GeometryStorage = FVGridGeometryStorage<
        CCTpfa::Detail::Face<GV>,
        CCTpfa::Detail::ScvWithId<GV, ThisType>,
        CCTpfa::Detail::ScvfWithId<GV, ThisType>,
        StorageSizes
    >;

    using FaceSeed = typename GG::FaceSeed;

    struct FaceNeighbors
    {
        LocalFaceNeighbors neighborScvs;
        LocalFaceNeighbors neighborScvfs;
    };

    struct ScvfLocalKey
    {
        std::size_t faceIndex;
        unsigned int idxInNeighbors;
    };

public:
    using GridGeometry = GG;
    using Element = typename GV::template Codim<0>::Entity;
    using SubControlVolume = typename GeometryStorage::SubControlVolume;
    using SubControlVolumeFace = typename GeometryStorage::SubControlVolumeFace;

    LocalView(const GridGeometry& gg)
    : gridGeometry_(gg)
    {}

    const GridGeometry& gridGeometry() const
    { return gridGeometry_; }

    void bind(const Element& e) &
    {
        clear_();
        set_(e);
        reserve_();
        makeInsideGeometries_();
        makeOutsideGeometries_();
    }

    void bindElement(const Element& e) &
    {
        clear_();
        set_(e);
        reserve_();
        makeInsideGeometries_();
    }

    LocalView bind(const Element& e) && { this->bind(e); return std::move(*this); }
    LocalView bindElement(const Element& e) && { this->bindElement(e); return std::move(*this); }

    friend std::ranges::range auto scvs(const LocalView& lv)
    { return std::views::single(lv.storage_.scv(0)); }

    friend std::ranges::range auto scvfs(const LocalView& lv)
    {
        return std::views::transform(
            std::views::iota(std::size_t{0}, lv.numInsideScvfs_),
            [&] (const std::size_t i) -> const SubControlVolumeFace& {
                return lv.storage_.scvf(i);
            }
        );
    }

    friend std::ranges::range auto neighborScvs(const LocalView& lv)
    {
        return std::views::transform(
            std::views::iota(std::size_t{1}, lv.storage_.numScvs()),
            [&] (const std::size_t i) -> const SubControlVolume& {
                return lv.storage_.scv(i);
            }
        );
    }

    friend std::ranges::range auto neighborScvfs(const LocalView& lv)
    {
        return std::views::transform(
            std::views::iota(lv.numInsideScvfs_, lv.storage_.numScvfs()),
            [&] (const std::size_t i) -> const SubControlVolumeFace& {
                return lv.storage_.scvf(i);
            }
        );
    }

    const Element& element() const
    { return *element_; }

    auto geometry(const SubControlVolume& scv) const
    { return gridGeometry_.element(scv.dofIndex()).geometry(); }

    auto geometry(const SubControlVolumeFace& scvf) const
    {
        const auto& facet = scvfToSeedMap_[scvf.id]->insideFacet();
        return gridGeometry_.element(facet.elementIndex).template subEntity<1>(facet.facetIndex).geometry();
    }

    bool onBoundary(const SubControlVolumeFace& scvf) const
    { return scvfToSeedMap_[scvf.id]->onBoundary(); }

    const SubControlVolume& insideScv(const SubControlVolumeFace& scvf) const
    {
        const auto [faceIdx, idxInNeighbors] = scvfLocalKeys_[scvf.id];
        return storage_.scv(faceNeighbors_[faceIdx].neighborScvs[idxInNeighbors]);
    }

    const SubControlVolume& outsideScv(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto [faceIdx, idxInNeighbors] = scvfLocalKeys_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < faceNeighbors_[faceIdx].neighborScvs.size());
        return storage_.scv(faceNeighbors_[faceIdx].neighborScvs[i]);
    }

    const SubControlVolumeFace& flipScvf(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto [faceIdx, idxInNeighbors] = scvfLocalKeys_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < faceNeighbors_[faceIdx].neighborScvfs.size());
        return storage_.scvf(faceNeighbors_[faceIdx].neighborScvfs[i]);
    }

private:
    void clear_()
    {
        element_.reset();
        storage_.clear();
        faceNeighbors_.clear();
        scvfLocalKeys_.clear();
        scvfToSeedMap_.clear();
    }

    void reserve_()
    {
        std::size_t numFaces = 0;
        for (unsigned int facetIdx = 0; facetIdx < element_->subEntities(1); ++facetIdx)
            for ([[maybe_unused]] const FaceSeed& faceSeed : gridGeometry_.faceSeeds_(eIdx_, facetIdx))
                numFaces++;

        std::size_t numScvs = 0;
        std::size_t numScvfs = 0;
        if constexpr (!Size::isDynamic<Traits::maxNumBranchesPerScvf>)
        {
            numScvs = maxElementStencilSize;
            numScvfs = numFaces*Traits::maxNumBranchesPerScvf;
        }
        else
        {
            numScvs = numFaces + 1;
            numScvfs = numFaces*2;
        }

        storage_.reserve({
            .numFaces = numFaces,
            .numScvfs = numScvfs,
            .numScvs = numScvfs + 1
        });
    }

    void set_(const Element& e)
    {
        element_ = e;
        eIdx_ = gridGeometry_.elementMapper().index(e);
    }

    void makeInsideGeometries_()
    {
        const auto& insideElemGeometry = element_->geometry();
        pushScv_(eIdx_, insideElemGeometry);
        assert(storage_.numScvs() == 1 && "storage has not been cleared beforehand");

        for (unsigned int facetIdx = 0; facetIdx < element_->subEntities(1); ++facetIdx)
        {
            for (const FaceSeed& faceSeed : gridGeometry_.faceSeeds_(eIdx_, facetIdx))
            {
                const auto& [elemIdx, elemFacetIdx] = faceSeed.insideFacet();
                const auto& element = gridGeometry_.element(elemIdx);
                const auto& faceGeo = element.template subEntity<1>(elemFacetIdx).geometry();
                const std::size_t faceIdx = storage_.pushFace({convexPolytopeVolume(faceGeo), faceGeo.center()});
                const std::size_t curScvfIdx = storage_.pushScvf({
                    storage_.face(faceIdx),
                    getNormal(insideElemGeometry, facetIdx)
                });
                storage_.scvf(curScvfIdx).id = curScvfIdx;

                static constexpr unsigned int insideScvIdx = 0;
                static constexpr unsigned int idxInFaceNeighbors = 0;
                faceNeighbors_.push_back({{insideScvIdx}, {static_cast<unsigned>(curScvfIdx)}});
                scvfLocalKeys_.push_back({faceIdx, idxInFaceNeighbors});
                scvfToSeedMap_.push_back(&faceSeed);
            }
        }

        numInsideScvfs_ = storage_.numScvfs();
    }

    void makeOutsideGeometries_()
    {
        for (unsigned int faceIdx = 0; faceIdx < numInsideScvfs_; ++faceIdx)
        {
            const auto& faceSeed = *scvfToSeedMap_[faceIdx];
            const auto& face = storage_.face(faceIdx);
            auto& faceNeighborScvs = faceNeighbors_[faceIdx].neighborScvs;
            auto& faceNeighborScvfs = faceNeighbors_[faceIdx].neighborScvfs;

            unsigned int faceNeighborIdx = 1;
            for (const auto& [outsideElemIdx, outsideFacetIdx] : faceSeed.facets())
            {
                if (outsideElemIdx == eIdx_)
                    continue;

                const auto& outsideElement = gridGeometry_.element(outsideElemIdx);
                const auto& outsideElementGeo = outsideElement.geometry();
                const auto outsideScvId = pushScv_(outsideElemIdx, outsideElementGeo);
                const auto outsideScvfId = storage_.pushScvf({
                    face,
                    getNormal(outsideElementGeo, outsideFacetIdx)
                });
                storage_.scvf(outsideScvfId).id = outsideScvfId;

                faceNeighborScvs.push_back(outsideScvId);
                faceNeighborScvfs.push_back(outsideScvfId);
                scvfLocalKeys_.push_back({faceIdx, faceNeighborIdx});
                scvfToSeedMap_.push_back(&faceSeed);
                faceNeighborIdx++;
            }
        }
    }

    auto pushScv_(const std::integral auto dofIndex, const typename Element::Geometry& geo)
    {
        const std::size_t scvIdx = storage_.pushScv({
            dofIndex,
            convexPolytopeVolume(geo),
            geo.center()
        });
        storage_.scv(scvIdx).id = scvIdx;
        return scvIdx;
    }

    const GridGeometry& gridGeometry_;

    std::optional<Element> element_;
    typename GV::IndexSet::IndexType eIdx_;
    std::size_t numInsideScvfs_;

    GeometryStorage storage_;
    DefaultStorage<FaceNeighbors, maxNumElementFaces> faceNeighbors_;
    DefaultStorage<ScvfLocalKey, maxNumStencilScvfs> scvfLocalKeys_;
    DefaultStorage<const FaceSeed*, maxNumStencilScvfs> scvfToSeedMap_;
};

} // end namespace Dumux

#endif
