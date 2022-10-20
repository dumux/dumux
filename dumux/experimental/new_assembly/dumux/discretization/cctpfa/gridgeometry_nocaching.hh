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

#include <utility>
#include <cstdint>
#include <cassert>
#include <optional>

#include <dune/common/float_cmp.hh>

#include <dumux/geometry/diameter.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/basegridgeometry.hh>

#include <dumux/experimental/new_assembly/dumux/common/storage.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/indexedintersections.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/fvgridgeometrystorage.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>

#include "detail.hh"

namespace Dumux {

//! Specialization for the case of caching geometries only locally
template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, false, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = CCTpfaGridGeometry<GV, false, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;

public:
    using GridView = GV;
    using Extrusion = Extrusion_t<Traits>;

    class LocalView;
    using SubControlVolume = typename LocalView::SubControlVolume;
    using SubControlVolumeFace = typename LocalView::SubControlVolumeFace;

    using DiscretizationMethod = DiscretizationMethods::CCTpfa;
    static constexpr DiscretizationMethod discMethod{};
    static constexpr auto maxElementStencilSize = CCTpfa::Detail::maxElementStencilSize<Traits>;

    explicit CCTpfaGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    { CCTpfa::Detail::checkOverlapSize(gridView); }

    //! Return the total number of degrees of freedom
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    friend LocalView localView(const CCTpfaGridGeometry& gg)
    { return {gg}; }
};

template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, false, Traits>::LocalView
{
    using GG = CCTpfaGridGeometry<GV, false, Traits>;
    using ThisType = typename GG::LocalView;

    static constexpr auto maxNumFaceNeighbors = CCTpfa::Detail::maxNumFaceNeighbors<Traits>;
    static constexpr auto maxNumStencilScvfs = CCTpfa::Detail::maxNumStencilScvfs<Traits>;
    static constexpr auto maxNumElementFaces = CCTpfa::Detail::maxNumElementFaces<Traits>;
    static constexpr auto maxElementStencilSize = CCTpfa::Detail::maxElementStencilSize<Traits>;
    static constexpr bool isNetworkGrid = int{GV::dimension} < int{GV::dimensionworld};

    struct StorageSizes
    {
        static constexpr auto maxNumScv = maxElementStencilSize;
        static constexpr auto maxNumScvf = maxNumStencilScvfs;
        static constexpr auto maxNumFaces = maxNumElementFaces;
    };

    using IndexedIntersections = Dumux::IndexedIntersections<GV, maxNumStencilScvfs>;
    using GeometryStorage = FVGridGeometryStorage<
        CCTpfa::Detail::Face<GV>,
        CCTpfa::Detail::ScvWithId<GV, ThisType>,
        CCTpfa::Detail::ScvfWithId<GV, ThisType>,
        StorageSizes
    >;

    struct LocalFaceConnectivity
    {
        bool boundary;
        DefaultStorage<std::size_t, maxNumFaceNeighbors> neighborScvs;
        DefaultStorage<std::size_t, maxNumFaceNeighbors> neighborScvfs;
    };

    struct ScvfLocalKey
    {
        std::size_t faceIndex;
        std::size_t idxInNeighbors;
    };

public:
    using GridGeometry = GG;
    using Element = typename GV::template Codim<0>::Entity;
    using SubControlVolume = typename GeometryStorage::SubControlVolume;
    using SubControlVolumeFace = typename GeometryStorage::SubControlVolumeFace;

    LocalView(const GridGeometry& gg)
    : gridGeometry_(gg)
    , indexedIntersections_(gg.gridView())
    {}

    const GridGeometry& gridGeometry() const
    { return gridGeometry_; }

    void bind(const Element& e) &
    {
        clear_();
        set_(e);
        reserve_();
        makeGeometries_<true>();
    }

    void bindElement(const Element& e) &
    {
        clear_();
        set_(e);
        reserve_();
        makeGeometries_<false>();
    }

    LocalView bind(const Element& e) && { this->bind(e); return std::move(*this); }
    LocalView bindElement(const Element& e) && { this->bindElement(e); return std::move(*this); }

    const Element& element() const
    { return *element_; }

    auto geometry(const SubControlVolume& scv) const
    { return gridGeometry_.element(scv.dofIndex()).geometry(); }

    auto geometry(const SubControlVolumeFace& scvf) const
    {
        const auto scvfFaceIdx = scvfLocalKeys_[scvf.id].faceIndex;
        for (const auto& [faceIdx, is] : IndexedIntersections{gridGeometry_.gridView(), *element_})
            if (faceIdx == scvfFaceIdx)
                return is.geometry();
        DUNE_THROW(Dune::InvalidStateException, "Could not get scvf geometry");
    }

    bool onBoundary(const SubControlVolumeFace& scvf) const
    {
        const auto scvfFaceIdx = scvfLocalKeys_[scvf.id].faceIndex;
        return faceConnectivity_[scvfFaceIdx].boundary;
    }

    std::size_t numOutsideNeighbors(const SubControlVolumeFace& scvf) const
    {
        const auto scvfFaceIdx = scvfLocalKeys_[scvf.id].faceIndex;
        return faceConnectivity_[scvfFaceIdx].neighborScvs.size() - 1;
    }

    const SubControlVolume& insideScv(const SubControlVolumeFace& scvf) const
    {
        const auto [faceIdx, idxInNeighbors] = scvfLocalKeys_[scvf.id];
        return storage_.scv(faceConnectivity_[faceIdx].neighborScvs[idxInNeighbors]);
    }

    const SubControlVolume& outsideScv(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto [faceIdx, idxInNeighbors] = scvfLocalKeys_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < faceConnectivity_[faceIdx].neighborScvs.size());
        return storage_.scv(faceConnectivity_[faceIdx].neighborScvs[i]);
    }

    const SubControlVolumeFace& flipScvf(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto [faceIdx, idxInNeighbors] = scvfLocalKeys_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < faceConnectivity_[faceIdx].neighborScvfs.size());
        return storage_.scvf(faceConnectivity_[faceIdx].neighborScvfs[i]);
    }

    friend std::ranges::range auto scvs(const LocalView& lv)
    { return std::views::single(lv.storage_.scv(0)); }

    friend std::ranges::range auto scvfs(const LocalView& lv)
    {
        return std::views::transform(
            lv.insideScvfIndices_,
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
            lv.neighborScvfIndices_,
            [&] (const std::size_t i) -> const SubControlVolumeFace& {
                return lv.storage_.scvf(i);
            }
        );
    }

private:
    template<bool withNeighbors>
    void update_(const Element& e)
    {
        clear_();
        set_(e);
        reserve_();
        makeGeometries_<withNeighbors>();
    }

    void clear_()
    {
        element_.reset();
        storage_.clear();
        insideScvfIndices_.clear();
        neighborScvfIndices_.clear();
        faceConnectivity_.clear();
        scvfLocalKeys_.clear();
    }

    void set_(const Element& e)
    {
        element_ = e;
    }

    void reserve_()
    {
        if constexpr (Size::isDynamic<maxNumStencilScvfs>)
        {
            // We reserve memory for speed-up and in order to avoid memory reallocation
            // since the scvfs store references to faces (reallocation would cause dangling refs)
            const auto isections = intersections(gridGeometry_.gridView(), *element_);
            const std::size_t numNeighbors = std::accumulate(
                isections.begin(),
                isections.end(),
                std::size_t{0},
                [] (const std::size_t current, const auto& is) {
                    return current + 1;
                }
            );
            storage_.reserve({
                .numFaces = numNeighbors,
                .numScvfs = numNeighbors*2 + 1,
                .numScvs = numNeighbors + 1
            });
        }
    }

    template<bool makeNeighborGeometries>
    void makeGeometries_()
    {
        assert(storage_.numScvs() == 0 && "Storage has not been cleared beforehand");
        assert(storage_.numScvfs() == 0 && "Storage has not been cleared beforehand");

        indexedIntersections_.set(*element_);
        pushScv_(gridGeometry_.elementMapper().index(*element_), element_->geometry());
        for (const auto& [fIdx, intersection] : indexedIntersections_)
        {
            const auto [faceIdx, inserted] = pushFace_(intersection, fIdx);
            if (inserted)
                pushInsideScvf_(faceIdx, intersection);

            if constexpr (makeNeighborGeometries)
            {
                if (intersection.neighbor())
                {
                    const auto& neighbor = intersection.outside();
                    const auto neighborIdx = gridGeometry_.elementMapper().index(neighbor);
                    const auto neighborScvIdx = pushScv_(neighborIdx, neighbor.geometry());
                    pushOutsideScvf_(faceIdx, neighborScvIdx, intersection);
                }
            }
        }
    }

    template<typename Intersection>
    std::pair<std::size_t, bool> pushFace_(const Intersection& is, int faceIdx) requires(!isNetworkGrid)
    {
        const auto& faceGeo = is.geometry();
        return {storage_.pushFace({faceGeo.volume(), faceGeo.center()}), true};
    }

    template<typename Intersection>
    std::pair<std::size_t, bool> pushFace_(const Intersection& is, int faceIdx) requires(isNetworkGrid)
    {
        if (faceExists_(faceIdx))
            return {faceIdx, false};
        const auto& faceGeo = is.geometry();
        return {storage_.pushFace({faceGeo.volume(), faceGeo.center()}), true};
    }

    bool faceExists_(int faceIdx) const
    { return faceIdx < storage_.numFaces(); }

    std::size_t pushScv_(const std::integral auto dofIndex, const typename Element::Geometry& geo)
    {
        const std::size_t scvIdx = storage_.pushScv({dofIndex, geo.volume(), geo.center()});
        storage_.scv(scvIdx).id = scvIdx;
        return scvIdx;
    }

    template<typename Face, typename GlobalCoordinate>
    std::size_t pushScvf_(const Face& face, GlobalCoordinate&& normal)
    {
        const std::size_t scvfIdx = storage_.pushScvf({face, std::forward<GlobalCoordinate>(normal)});
        storage_.scvf(scvfIdx).id = scvfIdx;
        return scvfIdx;
    }

    template<typename Intersection>
    std::size_t pushInsideScvf_(std::size_t faceIdx, const Intersection& is)
    {
        static constexpr unsigned int insideScvIdx = 0;
        static constexpr unsigned int idxInFaceNeighbors = 0;
        const std::size_t scvfIdx = pushScvf_(storage_.face(faceIdx), is.centerUnitOuterNormal());
        faceConnectivity_.push_back({is.boundary(), {insideScvIdx}, {scvfIdx}});
        insideScvfIndices_.push_back(scvfIdx);
        scvfLocalKeys_.push_back({faceIdx, idxInFaceNeighbors});
        return scvfIdx;
    }

    template<typename Intersection>
    std::size_t pushOutsideScvf_(std::size_t faceIdx,
                                 std::size_t insideScvIdx,
                                 const Intersection& intersection)
    {
        const std::size_t scvfIdx = pushScvf_(storage_.face(faceIdx), computeFlippedNormal_(intersection));
        const std::size_t idxInFaceNeighbors = faceConnectivity_[faceIdx].neighborScvs.size();
        faceConnectivity_[faceIdx].neighborScvs.push_back(insideScvIdx);
        faceConnectivity_[faceIdx].neighborScvfs.push_back(scvfIdx);
        neighborScvfIndices_.push_back(scvfIdx);
        scvfLocalKeys_.push_back({faceIdx, idxInFaceNeighbors});
        return scvfIdx;
    }

    template<typename Intersection>
    auto computeFlippedNormal_(const Intersection& is) requires(!isNetworkGrid)
    {
        auto n = is.centerUnitOuterNormal();
        n *= -1.0;
        return n;
    }

    template<typename Intersection>
    auto computeFlippedNormal_(const Intersection& is) requires(isNetworkGrid)
    {
        const auto isGeo = is.geometry();
        const auto isCenter = isGeo.center();
        const auto tolerance = 1e-6*diameter(isGeo);
        for (const auto& outsideIS : intersections(gridGeometry_.gridView(), is.outside()))
            if (Dune::FloatCmp::eq(isCenter, outsideIS.geometry().center(), tolerance))
                return outsideIS.centerUnitOuterNormal();
        DUNE_THROW(Dune::InvalidStateException, "Could not find the flip scvf normal");
    }

    const GridGeometry& gridGeometry_;
    IndexedIntersections indexedIntersections_;

    std::optional<Element> element_;
    GeometryStorage storage_;
    DefaultStorage<std::size_t, maxNumElementFaces> insideScvfIndices_;
    DefaultStorage<std::size_t, maxNumElementFaces> neighborScvfIndices_;
    DefaultStorage<LocalFaceConnectivity, maxNumElementFaces> faceConnectivity_;
    DefaultStorage<ScvfLocalKey, maxNumStencilScvfs> scvfLocalKeys_;
};

} // end namespace Dumux

#endif
