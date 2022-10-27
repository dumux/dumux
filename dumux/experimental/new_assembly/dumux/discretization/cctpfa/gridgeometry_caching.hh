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

#include <utility>
#include <cstdint>
#include <cassert>
#include <optional>
#include <limits>

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

//! Specialization for the case of caching geometries globally
template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, true, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = CCTpfaGridGeometry<GV, true, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;

    using Element = typename GV::template Codim<0>::Entity;
    using GeometryStorage = FVGridGeometryStorage<
        CCTpfa::Detail::Face<GV>,
        CCTpfa::Detail::ScvWithId<GV, ThisType>,
        CCTpfa::Detail::ScvfWithId<GV, ThisType>
    >;

    static constexpr auto unassignedIndex = std::numeric_limits<std::size_t>::max();
    static constexpr auto maxNumElementFaces = CCTpfa::Detail::maxNumElementFaces<Traits>;
    static constexpr auto maxNumFaceNeighbors = CCTpfa::Detail::maxNumFaceNeighbors<Traits>;
    static constexpr auto maxNumStencilScvfs = CCTpfa::Detail::maxNumStencilScvfs<Traits>;
    static constexpr bool isNetworkGrid = int{GV::dimension} < int{GV::dimensionworld};
    using IndexedIntersections = Dumux::IndexedIntersections<GV, maxNumStencilScvfs>;

    struct FaceConnectivity
    {
        bool boundary;
        DefaultStorage<std::size_t, maxNumFaceNeighbors> neighborScvs;
        DefaultStorage<std::size_t, maxNumFaceNeighbors> neighborScvfs;
    };

    struct ScvfKey
    {
        std::size_t faceIndex;
        std::uint_least8_t indexInNeighbors;
    };

public:
    using GridView = GV;
    using Extrusion = Extrusion_t<Traits>;

    class LocalView;
    using SubControlVolume = typename GeometryStorage::SubControlVolume;
    using SubControlVolumeFace = typename GeometryStorage::SubControlVolumeFace;

    using DiscretizationMethod = DiscretizationMethods::CCTpfa;
    static constexpr DiscretizationMethod discMethod{};
    static constexpr auto maxElementStencilSize = CCTpfa::Detail::maxElementStencilSize<Traits>;

    explicit CCTpfaGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    {
        CCTpfa::Detail::checkOverlapSize(gridView);
        update_();
    }

    //! Update all geometries (e.g. after grid adaption)
    void update(const GridView& gridView)
    {
        ParentType::update(gridView);
        update_();
    }

    //! Return the total number of degrees of freedom
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    friend LocalView localView(const CCTpfaGridGeometry& gg)
    { return {gg}; }

private:
    void update_()
    {
        clear_();
        reserve_();
        makeGeometries_();
        shrinkToFit_();

        assert(std::ranges::none_of(elementFaces_, [] (const auto& faceIndices) {
            return std::ranges::any_of(faceIndices, [] (const auto faceIndex) {
                return faceIndex == unassignedIndex;
            });
        }));
    }

    void clear_()
    {
        elementFaces_.clear();
        faceConnectivity_.clear();
        scvfKeys_.clear();
        storage_.clear();
    }

    void reserve_()
    {
        std::size_t numScvf = 0;
        std::size_t numBoundaryScvf = 0;

        IndexedIntersections indexedIntersections{this->gridView()};
        std::vector<int> elementFaceCounts(this->gridView().size(0), 0);
        std::ranges::for_each(elements(this->gridView()), [&] (const auto& element) {
            indexedIntersections.set(element);
            const auto eIdx = this->elementMapper().index(element);
            for (const auto& [fIdx, is] : indexedIntersections)
            {
                elementFaceCounts[eIdx] = std::max(elementFaceCounts[eIdx], fIdx + 1);
                numScvf++;
                if (is.boundary())
                    numBoundaryScvf++;
            };
        });

        elementFaces_.resize(this->gridView().size(0));
        for (std::size_t eIdx = 0; eIdx < elementFaceCounts.size(); ++eIdx)
            elementFaces_[eIdx].resize(elementFaceCounts[eIdx], unassignedIndex);

        const std::size_t numFaces = std::size_t{(numScvf - numBoundaryScvf)/2} + numBoundaryScvf;
        faceConnectivity_.reserve(numFaces);
        scvfKeys_.reserve(numScvf);
        storage_.reserve({
            .numFaces = numFaces,
            .numScvfs = numScvf,
            .numScvs = static_cast<std::size_t>(this->gridView().size(0))
        });
    }

    void makeGeometries_()
    {
        IndexedIntersections indexedIntersections{this->gridView()};
        for (typename GV::IndexSet::IndexType eIdx = 0; eIdx < this->gridView().size(0); ++eIdx)
        {
            const auto& element = this->element(eIdx);
            const auto& elemGeom = element.geometry();
            const auto tolerance = 1e-6*diameter(elemGeom);
            pushScv_(eIdx, elemGeom);
            indexedIntersections.set(element);
            makeFacesAndScvfs_(indexedIntersections, eIdx, tolerance);
        }
    }

    void makeFacesAndScvfs_(const IndexedIntersections& indexedIntersections,
                            std::integral auto eIdx,
                            std::floating_point auto tolerance) requires(!isNetworkGrid)
    {
        for (const auto& [localFaceIdx, is] : indexedIntersections)
        {
            const auto& isGeo = is.geometry();
            if (const auto faceIdx = getFaceIndex_(is, isGeo, tolerance); faceIdx)
            {
                elementFaces_[eIdx][localFaceIdx] = *faceIdx;
                pushScvf_(*faceIdx, eIdx, is.centerUnitOuterNormal());
            }
            else
                pushScvf_(
                    *pushFace_(is, isGeo, eIdx, localFaceIdx),
                    eIdx,
                    is.centerUnitOuterNormal()
                );
        }
    }

    void makeFacesAndScvfs_(const IndexedIntersections& indexedIntersections,
                            std::integral auto eIdx,
                            std::floating_point auto tolerance) requires(isNetworkGrid)
    {
        for (const auto& [localFaceIdx, is] : indexedIntersections)
        {
            const auto& isGeo = is.geometry();
            if (const auto faceIdx = getFaceIndex_(is, isGeo, tolerance);
                faceIdx && elementFaces_[eIdx][localFaceIdx] == unassignedIndex)
            {
                elementFaces_[eIdx][localFaceIdx] = *faceIdx;
                pushScvf_(*faceIdx, eIdx, is.centerUnitOuterNormal());
            }
            else if (const auto newFaceIdx = pushFace_(is, isGeo, eIdx, localFaceIdx); newFaceIdx)
                pushScvf_(*newFaceIdx, eIdx, is.centerUnitOuterNormal());
        }
    }

    void shrinkToFit_()
    {
        storage_.shrinkToFit();
        faceConnectivity_.shrink_to_fit();
        scvfKeys_.shrink_to_fit();
    }

    template<typename Intersection>
    std::optional<std::size_t> getFaceIndex_(const Intersection& is,
                                             const typename Intersection::Geometry& isGeo,
                                             const typename Intersection::ctype tolerance) const
    {
        if (!is.neighbor())
            return std::nullopt;
        const auto isCenter = isGeo.center();
        const auto isAssigned = [] (const std::size_t idx) { return idx != unassignedIndex; };
        const auto& outsideFaceIndices = elementFaces_[this->elementMapper().index(is.outside())];
        for (const auto fIdx : std::views::filter(outsideFaceIndices, isAssigned))
            if (Dune::FloatCmp::eq(storage_.face(fIdx).center(), isCenter, tolerance))
                return fIdx;
        return std::nullopt;
    }

    template<typename Intersection>
    std::optional<std::size_t> pushFace_(const Intersection& is,
                                         const typename Intersection::Geometry& isGeo,
                                         const std::integral auto eIdx,
                                         [[maybe_unused]] const int localFaceIdx) requires(!isNetworkGrid)
    {
        const std::size_t faceIdx = storage_.pushFace({isGeo.volume(), isGeo.center()});
        elementFaces_[eIdx][localFaceIdx] = faceIdx;
        faceConnectivity_.push_back({is.boundary(), {}, {}});
        return faceIdx;
    }

    template<typename Intersection>
    std::optional<std::size_t> pushFace_(const Intersection& is,
                                         const typename Intersection::Geometry& isGeo,
                                         const std::integral auto eIdx,
                                         const int localFaceIdx) requires(isNetworkGrid)
    {
        if (elementFaces_[eIdx][localFaceIdx] != unassignedIndex)
            return std::nullopt;
        const std::size_t faceIdx = storage_.pushFace({isGeo.volume(), isGeo.center()});
        elementFaces_[eIdx][localFaceIdx] = faceIdx;
        faceConnectivity_.push_back({is.boundary(), {}, {}});
        return faceIdx;
    }

    std::size_t pushScv_(const std::integral auto dofIndex, const typename Element::Geometry& geo)
    {
        const std::size_t scvIdx = storage_.pushScv({dofIndex, geo.volume(), geo.center()});
        storage_.scv(scvIdx).id = scvIdx;
        return scvIdx;
    }

    template<typename GlobalCoordinate>
    std::size_t pushScvf_(const std::integral auto faceIdx,
                          const std::integral auto insideScvIdx,
                          GlobalCoordinate&& normal)
    {
        const std::size_t scvfIdx = storage_.pushScvf({
            storage_.face(faceIdx),
            std::forward<GlobalCoordinate>(normal)
        });
        storage_.scvf(scvfIdx).id = scvfIdx;
        scvfKeys_.push_back({
            .faceIndex = faceIdx,
            .indexInNeighbors = static_cast<std::uint_least8_t>(
                faceConnectivity_[faceIdx].neighborScvs.size()
            )
        });
        faceConnectivity_[faceIdx].neighborScvs.push_back(insideScvIdx);
        faceConnectivity_[faceIdx].neighborScvfs.push_back(scvfIdx);
        return scvfIdx;
    }

    friend class LocalView;
    GeometryStorage storage_;
    std::vector<ScvfKey> scvfKeys_;
    std::vector<FaceConnectivity> faceConnectivity_;
    std::vector<DefaultStorage<std::size_t, maxNumElementFaces>> elementFaces_;
};

template<typename GV, typename Traits>
class CCTpfaGridGeometry<GV, true, Traits>::LocalView
{
    using GG = CCTpfaGridGeometry<GV, true, Traits>;
    using ThisType = typename GG::LocalView;

    static constexpr auto maxNumElementFaces = CCTpfa::Detail::maxNumElementFaces<Traits>;
    static constexpr auto maxNumStencilScvfs = CCTpfa::Detail::maxNumStencilScvfs<Traits>;

public:
    using GridGeometry = GG;
    using Element = typename GV::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    LocalView(const GridGeometry& gg)
    : gridGeometry_(gg)
    {}

    const GridGeometry& gridGeometry() const
    { return gridGeometry_; }

    void bind(const Element& e) &
    { update_<true>(e); }

    void bindElement(const Element& e) &
    { update_<false>(e); }

    LocalView bind(const Element& e) && { this->bind(e); return std::move(*this); }
    LocalView bindElement(const Element& e) && { this->bindElement(e); return std::move(*this); }

    const Element& element() const
    { return *element_; }

    auto geometry(const SubControlVolume& scv) const
    { return gridGeometry_.element(scv.dofIndex()).geometry(); }

    auto geometry(const SubControlVolumeFace& scvf) const
    {
        const auto [faceIndex, indexInNeighbors] = gridGeometry_.scvfKeys_[scvf.id];
        const auto eIdx = gridGeometry_.faceConnectivity_[faceIndex].neighborScvs[indexInNeighbors];
        const auto element = gridGeometry_.element(eIdx);
        for (const auto& [localFaceIdx, is] : IndexedIntersections{gridGeometry_.gridView(), element})
            if (gridGeometry_.elementFaces_[eIdx][localFaceIdx] == faceIndex)
                return is.geometry();
        DUNE_THROW(Dune::InvalidStateException, "Could not get scvf geometry");
    }

    bool onBoundary(const SubControlVolumeFace& scvf) const
    {
        const auto [faceIndex, _] = gridGeometry_.scvfKeys_[scvf.id];
        return gridGeometry_.faceConnectivity_[faceIndex].boundary;
    }

    std::size_t numOutsideNeighbors(const SubControlVolumeFace& scvf) const
    {
        const auto [faceIndex, _] = gridGeometry_.scvfKeys_[scvf.id];
        return gridGeometry_.faceConnectivity_[faceIndex].neighborScvs.size() - 1;
    }

    std::size_t indexInElement(const SubControlVolume& scv) const
    {
        assert(scv.dofIndex() == eIdx_ && "Index in element only available for scvs inside the bound element");
        return 0;
    }

    std::size_t indexInElement(const SubControlVolumeFace& scvf) const
    {
        const auto [faceIndex, _] = gridGeometry_.scvfKeys_[scvf.id];
        const auto& faces = gridGeometry_.elementFaces_[eIdx_];
        auto it = std::ranges::find(faces, faceIndex);
        assert(it != std::ranges::end(faces) && "Index in element only available for scvfs inside the bound element");
        return std::ranges::distance(std::ranges::begin(faces), it);
    }

    const SubControlVolume& insideScv(const SubControlVolumeFace& scvf) const
    {
        const auto [faceIndex, indexInNeighbors] = gridGeometry_.scvfKeys_[scvf.id];
        const auto eIdx = gridGeometry_.faceConnectivity_[faceIndex].neighborScvs[indexInNeighbors];
        return getScv_(eIdx);
    }

    const SubControlVolume& outsideScv(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto [faceIdx, idxInNeighbors] = gridGeometry_.scvfKeys_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < gridGeometry_.faceConnectivity_[faceIdx].neighborScvs.size());
        return getScv_(gridGeometry_.faceConnectivity_[faceIdx].neighborScvs[i]);
    }

    const SubControlVolumeFace& flipScvf(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    {
        const auto [faceIdx, idxInNeighbors] = gridGeometry_.scvfKeys_[scvf.id];
        if (i >= idxInNeighbors) i++;
        assert(i < gridGeometry_.faceConnectivity_[faceIdx].neighborScvfs.size());
        return getScvf_(gridGeometry_.faceConnectivity_[faceIdx].neighborScvfs[i]);
    }

    friend std::ranges::range auto scvs(const LocalView& lv)
    { return std::views::single(lv.getScv_(lv.eIdx_)); }

    friend std::ranges::range auto scvfs(const LocalView& lv)
    {
        return std::views::transform(
            lv.insideScvfIndices_,
            [&] (const std::size_t i) -> const SubControlVolumeFace& {
                return lv.getScvf_(i);
            }
        );
    }

    friend std::ranges::range auto neighborScvs(const LocalView& lv)
    {
        return std::views::transform(
            lv.neighborScvIndices_,
            [&] (const std::size_t i) -> const SubControlVolume& {
                return lv.getScv_(i);
            }
        );
    }

    friend std::ranges::range auto neighborScvfs(const LocalView& lv)
    {
        return std::views::transform(
            lv.neighborScvfIndices_,
            [&] (const std::size_t i) -> const SubControlVolumeFace& {
                return lv.getScvf_(i);
            }
        );
    }

private:
    const SubControlVolume& getScv_(std::size_t i) const
    { return gridGeometry_.storage_.scv(i); }

    const SubControlVolumeFace& getScvf_(std::size_t i) const
    { return gridGeometry_.storage_.scvf(i); }

    template<bool withNeighbors>
    void update_(const Element& e)
    {
        clear_();
        set_(e);
        reserve_<withNeighbors>();
        prepareScvfIndexSets_<withNeighbors>();
    }

    void clear_()
    {
        insideScvfIndices_.clear();
        neighborScvfIndices_.clear();
        neighborScvIndices_.clear();
    }

    void set_(const Element& e)
    {
        element_ = e;
        eIdx_ = gridGeometry_.elementMapper().index(e);
    }

    template<bool withNeighbors>
    void reserve_()
    {
        if constexpr (Size::isDynamic<maxNumElementFaces>)
            insideScvfIndices_.reserve(gridGeometry_.elementFaces_[eIdx_].size());

        if constexpr (withNeighbors && Size::isDynamic<maxNumStencilScvfs>)
        {
            const std::size_t numNeighbors = std::accumulate(
                gridGeometry_.elementFaces_[eIdx_].begin(),
                gridGeometry_.elementFaces_[eIdx_].end(),
                std::size_t{0},
                [&] (const std::size_t current, const std::size_t faceIdx) {
                    return current + gridGeometry_.faceConnectivity_[faceIdx].neighborScvs.size() - 1;
            });
            neighborScvIndices_.reserve(numNeighbors);
            neighborScvfIndices_.reserve(numNeighbors);
        }
    }

    template<bool withNeighbors>
    void prepareScvfIndexSets_()
    {
        std::ranges::for_each(gridGeometry_.elementFaces_[eIdx_], [&] (const auto& fIdx) {
            const auto& connectivity = gridGeometry_.faceConnectivity_[fIdx];
            for (std::size_t i = 0; i < connectivity.neighborScvs.size(); ++i)
            {
                if (connectivity.neighborScvs[i] == eIdx_)
                    insideScvfIndices_.push_back(connectivity.neighborScvfs[i]);
                else if constexpr (withNeighbors)
                {
                    neighborScvfIndices_.push_back(connectivity.neighborScvfs[i]);
                    neighborScvIndices_.push_back(connectivity.neighborScvs[i]);
                }
            }
        });
    }

    const GridGeometry& gridGeometry_;
    std::optional<Element> element_;
    typename GV::IndexSet::IndexType eIdx_;

    DefaultStorage<std::size_t, maxNumElementFaces> insideScvfIndices_;
    DefaultStorage<std::size_t, maxNumStencilScvfs> neighborScvfIndices_;
    DefaultStorage<std::size_t, maxNumStencilScvfs> neighborScvIndices_;
};

} // end namespace Dumux

#endif
