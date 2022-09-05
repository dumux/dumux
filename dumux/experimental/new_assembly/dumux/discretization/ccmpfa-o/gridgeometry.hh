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
 * \ingroup CCMpfaDiscretization
 * \brief TODO: Doc me
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_O_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCMPFA_O_GRID_GEOMETRY_HH

#include <vector>
#include <utility>
#include <concepts>
#include <algorithm>
#include <type_traits>
#include <unordered_map>
#include <optional>
#include <numeric>
#include <limits>

#include <dune/common/reservedvector.hh>

#include <dumux/common/defaultmappertraits.hh>
#include <dumux/geometry/volume.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/cellcentered/connectivitymap.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>
#include <dumux/discretization/extrusion.hh>


#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>

namespace Dumux {

// TODO: network grids
// TODO: flip scvf
// TODO: check overlap size

// TODO: proper helper (with compile-time overloads)
template<typename Geometry>
auto getNormal(const Geometry& geo, unsigned int facetIdx)
{
    static_assert(Geometry::mydimension == 2);
    const auto rotateClockWiseAndScale = [] (auto&& v) {
        std::swap(v[0], v[1]);
        v[1] *= -1.0;
        v /= v.two_norm();
        return v;
    };

    using Dune::referenceElement;
    const auto refElement = referenceElement(geo);
    if (geo.type() == Dune::GeometryTypes::quadrilateral)
    {
        switch (facetIdx)
        {
            case 0: return rotateClockWiseAndScale(
                geo.global(refElement.position(0, 2)) -
                geo.global(refElement.position(2, 2))
            );
            case 1: return rotateClockWiseAndScale(
                geo.global(refElement.position(3, 2)) -
                geo.global(refElement.position(1, 2))
            );
            case 2: return rotateClockWiseAndScale(
                geo.global(refElement.position(1, 2)) -
                geo.global(refElement.position(0, 2))
            );
            case 3: return rotateClockWiseAndScale(
                geo.global(refElement.position(2, 2)) -
                geo.global(refElement.position(3, 2))
            );
        }
    }

    DUNE_THROW(Dune::NotImplemented, "TODO");
}

/*!
 * \ingroup CCMpfaDiscretization
 * \todo Doc me
 */
template<typename GV,
         bool cache = true,
         typename Traits = DefaultMapperTraits<GV>>
class CCMpfaOGridGeometry
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = CCMpfaOGridGeometry<GV, cache, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;

    using LocalIndex = std::uint_least8_t;
    using GridIndex = typename GV::IndexSet::IndexType;
    using Element = typename GV::template Codim<0>::Entity;
    using Coordinate = typename Element::Geometry::GlobalCoordinate;
    using ctype = typename GV::ctype;

    static constexpr int dim = GV::dimension;
    static constexpr int dimWorld = GV::dimensionworld;
    static constexpr bool isNetworkGrid = dim < dimWorld;
    static constexpr GridIndex undefined = std::numeric_limits<GridIndex>::max();

    template<typename T>
    using ScvfNeighborStorage = std::conditional_t<isNetworkGrid,
                                                   std::vector<T>,
                                                   Dune::ReservedVector<T, 2>>;

    struct ScvfSeed
    {
        GridIndex elementIndex;
        LocalIndex facetIndex;
    };

    using ScvfsPerGridFace = ScvfNeighborStorage<ScvfSeed>;

    template<int codim, typename ID = LocalIndex>
    class SubControlEntity
    {
    public:
        SubControlEntity() = default;
        SubControlEntity(ID id) noexcept : id_(std::move(id)) {}
        const ID& id() const { return id_; }
    private:
        ID id_;
    };

public:
    using SubControlVolume = SubControlEntity<0>;
    using SubControlVolumeFace = SubControlEntity<1>;

    class LocalView
    {
        static constexpr int maxNumFacesPerCell = 8;  // cubes
        static constexpr int maxNumScvs = maxNumFacesPerCell + 1;
        static constexpr int maxNumScvfs = maxNumFacesPerCell*2;

        struct Face
        {
            ctype area;
            Coordinate center;
            ScvfNeighborStorage<LocalIndex> localNeighborScvIndices;
        };

        struct ScvStorage
        {
            SubControlVolume scv;
            GridIndex dofIndex;
            Coordinate center;
            ctype volume;
        };

        struct ScvfStorage
        {
            SubControlVolumeFace scvf;
            Coordinate normal;
            LocalIndex facetIdx;
            LocalIndex indexInFaceNeighbors;
        };

    public:
        using GridGeometry = ThisType;

        explicit LocalView(const GridGeometry& gg)
        : gridGeom_(gg)
        {}

        void bind(const Element& element) &
        { bind_(element); }

        LocalView bind(const Element& element) &&
        {
            element_ = element;
            return std::move(*this);
        }

        GridIndex dofIndex(const SubControlVolume& scv) const
        { return scvs_[scv.id()].dofIndex; }

        auto geometry(const SubControlVolume& scv) const
        {
            if (scv.id() == 0)
                return element_->geometry();
            return gridGeom_.element(scvs_[scv.id()].dofIndex).geometry();
        }

        auto geometry(const SubControlVolumeFace& scvf) const
        {
            const auto facetIdx = scvfs_[scvf.id()].facetIdx;
            return element_->template subEntity<1>(facetIdx).geometry();
        }

        ctype volume(const SubControlVolume& scv) const
        { return scvs_[scv.id()].volume; }

        ctype area(const SubControlVolumeFace& scvf) const
        { return getFace_(scvf).area; }

        const Coordinate& center(const SubControlVolume& scv) const
        { return scvs_[scv.id()].center; }

        const Coordinate& center(const SubControlVolumeFace& scvf) const
        { return getFace_(scvf).center; }

        const SubControlVolume& insideScv(const SubControlVolumeFace& scvf) const
        {
            const auto indexInNeighbors = scvfs_[scvf.id()].indexInFaceNeighbors;
            return scvs_[getFace_(scvf).localNeighborScvIndices[indexInNeighbors]].scv;
        }

        const SubControlVolume& outsideScv(const SubControlVolumeFace& scvf, unsigned int i) const
        {
            if (i >= scvfs_[scvf.id()].indexInFaceNeighbors)
                i++;
            return scvs_[getFace_(scvf).localNeighborScvIndices[i]].scv;
        }

        std::size_t numOutsideScvs(const SubControlVolumeFace& scvf) const
        { return getFace_(scvf).localNeighborScvIndices.size() - 1; }

        bool onBoundary(const SubControlVolumeFace& scvf) const
        { return numOutsideScvs(scvf) == 0; }

        const Coordinate& unitOuterNormal(const SubControlVolumeFace& scvf) const
        { return scvfs_[scvf.id()].normal; }

        friend auto scvs(const LocalView& lv)
        {
            return std::views::transform(std::views::iota(0, 1), [&] (int scvIdx) {
                return lv.scvs_[scvIdx].scv;
            });
        }

        friend auto scvfs(const LocalView& lv)
        {
            const int numFacets = static_cast<int>(lv.faces_.size());
            return std::views::transform(std::views::iota(0, numFacets), [&] (int fIdx) {
                return lv.scvfs_[fIdx].scvf;
            });
        }

    private:
        template<std::integral I>
        LocalIndex makeLocalIndex_(const I& i) const
        { return static_cast<LocalIndex>(i); }

        const Face& getFace_(const SubControlVolumeFace& scvf) const
        { return faces_[scvfs_[scvf.id()].facetIdx]; }

        void clear_()
        {
            scvs_.clear();
            scvfs_.clear();
            faces_.clear();
        }

        void bind_(const Element& element)
        {
            clear_();
            element_ = element;
            const auto& elemGeo = element.geometry();
            const auto insideElemIdx = gridGeom_.elementMapper().index(element);
            const auto insideScvLocalIdx = addScv_(elemGeo, insideElemIdx);
            makeFaces_(element, elemGeo, insideScvLocalIdx);
            makeOutsideEntities_();
        }

        void makeFaces_(const Element& element,
                        const typename Element::Geometry& geo,
                        LocalIndex insideScvLocalIdx)
        {
            const auto numFaces = static_cast<int>(element.subEntities(1));
            std::ranges::for_each(std::views::iota(0, numFaces), [&] (auto facetIdx) {
                const auto& faceGeom = element.template subEntity<1>(facetIdx).geometry();
                faces_.push_back(Face{
                    convexPolytopeVolume(faceGeom),
                    faceGeom.center(),
                    {insideScvLocalIdx}
                });
                scvfs_.push_back(ScvfStorage{
                    SubControlVolumeFace{makeLocalIndex_(scvfs_.size())},
                    getNormal(geo, facetIdx),
                    makeLocalIndex_(facetIdx),
                    LocalIndex{0}
                });
            });
        }

        void makeOutsideEntities_()
        {
            const auto insideElementIndex = scvs_[0].dofIndex;
            const int numFacets = static_cast<int>(faces_.size());
            std::ranges::for_each(std::views::iota(0, numFacets), [&] (auto localFaceId) {
                const auto faceIdx = gridGeom_.faceIndicesOfScv_[insideElementIndex][localFaceId];
                const auto& faceNeighborScvfs = gridGeom_.scvfsPerGridFace_[faceIdx];
                std::ranges::for_each(faceNeighborScvfs, [&] (const ScvfSeed& seed) {
                    if (seed.elementIndex != insideElementIndex)
                        addOutsideScvAndScvf_(seed.elementIndex,
                                              seed.facetIndex,
                                              makeLocalIndex_(localFaceId));
                });
            });
        }

        void addOutsideScvAndScvf_(GridIndex eIdx,
                                   LocalIndex facetIdx,
                                   LocalIndex localFaceId)
        {
            const auto element = gridGeom_.element(eIdx);
            const auto& elemGeo = element.geometry();
            const auto localScvIdx = addScv_(elemGeo, eIdx);

            auto& faceNeighbors = faces_[localFaceId].localNeighborScvIndices;
            const auto idxInNeighbors = faceNeighbors.size();
            faceNeighbors.push_back(localScvIdx);
            scvfs_.push_back(ScvfStorage{
                SubControlVolumeFace{makeLocalIndex_(scvfs_.size())},
                getNormal(elemGeo, facetIdx),
                makeLocalIndex_(localFaceId),
                makeLocalIndex_(idxInNeighbors)
            });
        }

        LocalIndex addScv_(const typename Element::Geometry& geo, GridIndex dofIndex)
        {
            const auto localIndex = static_cast<LocalIndex>(scvs_.size());
            scvs_.push_back(ScvStorage{
                SubControlVolume{localIndex},
                dofIndex,
                geo.center(),
                convexPolytopeVolume(geo)
            });
            return localIndex;
        }

        const GridGeometry& gridGeom_;
        std::optional<Element> element_;
        Dune::ReservedVector<ScvStorage, maxNumScvs> scvs_;
        Dune::ReservedVector<ScvfStorage, maxNumScvfs> scvfs_;
        Dune::ReservedVector<Face, maxNumFacesPerCell> faces_;
    };

    // using LocalView = CCMpfaOGridGeometryLocalView<ThisType>;
    // using SubControlVolume = CCMpfaOSubControlVolumeEntity<GridIndex, Coordinate, 0>;
    // using SubControlVolumeFace = ScvfView;
    using DofMapper = typename Traits::ElementMapper;
    // using Extrusion = Extrusion_t<Traits>;

    // //! export the discretization method this geometry belongs to
    // using DiscretizationMethod = DiscretizationMethods::CCTpfa;
    // static constexpr DiscretizationMethod discMethod{};

    // //! The maximum admissible stencil size (used for static memory allocation during assembly)
    // static constexpr int maxElementStencilSize = LocalView::maxNumElementScvfs*Traits::maxNumScvfNeighbors + 1;

    //! export the grid view type
    using GridView = GV;

    //! Constructor
    CCMpfaOGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    {
        // // Check if the overlap size is what we expect
        // if (!CheckOverlapSize<DiscretizationMethod>::isValid(gridView))
        //     DUNE_THROW(Dune::InvalidStateException, "The cctpfa discretization method needs at least an overlap of 1 for parallel computations. "
        //                                              << " Set the parameter \"Grid.Overlap\" in the input file.");

        update_();
    }

    friend auto localView(const ThisType& gg)
    { return LocalView{gg}; }

    //! Here, the dof mapper is the element mapper.
    //! This is to have a better chance to have the same main files for different fv schemes
    const DofMapper& dofMapper() const
    { return this->elementMapper(); }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    //! Return the total number of sub control volumes
    std::size_t numScv() const
    { return faceIndicesOfScv_.size(); }

    //! Return the total number of sub control volume faces
    std::size_t numScvf() const
    {
        return std::accumulate(
            scvfsPerGridFace_.begin(),
            scvfsPerGridFace_.end(),
            std::size_t{0},
            [] (std::size_t current, const auto& neighbors) {
                return current + neighbors.size();
        });
    }

    // //! Return the geometry of an scv
    // auto geometry(const SubControlVolume& scv) const
    // { return this->element(scv.id()).geometry(); }

    // //! Return the dof index associated with an scv
    // GridIndex dofIndex(const SubControlVolume& scv) const
    // { return scv.id(); }

    // //! Return the geometry of an scvf
    // auto geometry(const SubControlVolumeFace& scvf) const
    // {
    //     const auto [index, eIdx, facetIdx] = scvf.id();
    //     return this->element(eIdx).template subEntity<1>(facetIdx).geometry();
    // }

    // //! Return true if the given scvf is on the boundary
    // bool onBoundary(const SubControlVolumeFace& scvf) const
    // { return numOutsideScvs(scvf) == 0; }

    // //! Return the number of "outside" scvs connected to an scvf
    // std::size_t numOutsideScvs(const SubControlVolumeFace& scvf) const
    // { return scvfNeighborIndices_[scvf.id().index].size() - 1; }

    // //! Return the "inside" scv of an scvf
    // const SubControlVolume& insideScv(const SubControlVolumeFace& scvf) const
    // {
    //     return scvs_[
    //         scvfNeighborIndices_[scvf.id().index]
    //                             [scvf.indexInNeighbors()]
    //     ];
    // }

    // //! Return the i-th "outside" scv of an scvf
    // const SubControlVolume& outsideScv(const SubControlVolumeFace& scvf, unsigned int i = 0) const
    // {
    //     assert(i < scvfNeighborIndices_[scvf.id().index].size() - 1);
    //     if (i >= scvf.indexInNeighbors())
    //         i++;
    //     return scvs_[scvfNeighborIndices_[scvf.id().index][i]];
    // }

    // //! Return the normal vector on an scvf
    // Coordinate unitOuterNormal(const SubControlVolumeFace& scvf) const
    // {
    //     static_assert(!isNetworkGrid);
    //     const auto [index, eIdx, facetIdx] = scvf.id();
    //     const auto& geo = this->element(eIdx).geometry();
    //     auto n = getNormal(geo, facetIdx);
    //     if (scvf.isOutsideView())
    //         n *= -1.0;
    //     return n;
    // }

    // //! Iterator range over all scvs
    // friend const std::vector<SubControlVolume>& scvs(const ThisType& gg) requires(cache)
    // { return gg.scvs_; }

    // //! Iterator range over the scvfs of an scv
    // friend auto scvfs(const ThisType& gg, const SubControlVolume& scv) requires(cache)
    // {
    //     return std::views::transform(
    //         gg.scvfIndicesOfScv_[gg.dofIndex(scv)],
    //         [&] (const auto& scvScvfId) {
    //             return ScvfView{
    //                 gg.scvfs_[scvScvfId.scvfIndex],
    //                 scvScvfId.myIndexInScvfNeighbors
    //             };
    //         }
    //     );
    // }

    // //! Iterator range over all scvs of an element
    // friend auto scvs(const ThisType& gg, const Element& element) requires(cache)
    // {
    //     const auto eIdx = gg.elementMapper().index(element);
    //     return std::views::transform(
    //         std::views::iota(eIdx, eIdx+1),
    //         [&] (auto eIdx) -> const SubControlVolume& {
    //             return gg.scvs_[eIdx];
    //         }
    //     );
    // }

    // //! Iterator range over all scvfs of an element
    // friend auto scvfs(const ThisType& gg, const Element& element) requires(cache)
    // {
    //     const auto eIdx = gg.elementMapper().index(element);
    //     const auto& scv = gg.scvs_[eIdx];
    //     return scvfs(gg, scv);
    // }

private:
    friend class LocalView;

    void reset_()
    {
        faceIndicesOfScv_.clear();
        faceIndicesOfScv_.resize(this->gridView().size(0));
        for (const auto& e : elements(this->gridView()))
            faceIndicesOfScv_[this->elementMapper().index(e)].resize(
                e.subEntities(1), undefined
            );

        scvfsPerGridFace_.clear();
        scvfsPerGridFace_.reserve(this->gridView().size(1));
    }

    void update_()
    {
        reset_();
        for (const auto& element : elements(this->gridView()))
            for (const auto& is : intersections(this->gridView(), element))
                registerFace_(is, this->elementMapper().index(element));

        assert(std::ranges::all_of(faceIndicesOfScv_, [] (const auto& indices) {
            return std::ranges::all_of(indices, [] (auto idx) {
                return idx != undefined;
            });
        }));
    }

    void registerFace_(const typename GV::Intersection& is, const GridIndex insideElemIdx)
    {
        if (!is.neighbor())
        {
            scvfsPerGridFace_.push_back({});
            addNeighbor_(scvfsPerGridFace_.size() - 1, insideElemIdx, is.indexInInside());
        }
        else
        {
            const auto outside = is.outside();
            const auto outsideElemIdx = this->elementMapper().index(outside);
            const auto insideFaceIdx = getFaceIndex_(insideElemIdx, is.indexInInside());
            const auto outsideFaceIdx = getFaceIndex_(outsideElemIdx, is.indexInOutside());
            assert((!insideFaceIdx || !outsideFaceIdx) || *insideFaceIdx == *outsideFaceIdx);

            if (insideFaceIdx && outsideFaceIdx)
                return;

            if (insideFaceIdx && !outsideFaceIdx)
                addNeighbor_(*insideFaceIdx, outsideElemIdx, is.indexInOutside());
            else if (outsideFaceIdx && !insideFaceIdx)
                addNeighbor_(*outsideFaceIdx, insideElemIdx, is.indexInInside());
            else
            {
                scvfsPerGridFace_.push_back({});
                addNeighbor_(scvfsPerGridFace_.size() - 1, insideElemIdx, is.indexInInside());
                addNeighbor_(scvfsPerGridFace_.size() - 1, outsideElemIdx, is.indexInOutside());
            }
        }
    }

    void addNeighbor_(const GridIndex faceIdx,
                      const GridIndex elementIndex,
                      const LocalIndex facetIndex)
    {
        scvfsPerGridFace_[faceIdx].emplace_back(elementIndex, facetIndex);
        faceIndicesOfScv_[elementIndex][facetIndex] = faceIdx;
    }

    std::optional<GridIndex> getFaceIndex_(const GridIndex elemIdx,
                                           const unsigned int facetIdx) const
    {
        if (!faceIndicesOfScv_[elemIdx].empty())
            if (faceIndicesOfScv_[elemIdx][facetIdx] != undefined)
                return faceIndicesOfScv_[elemIdx][facetIdx];
        return {};
    }

    std::vector<ScvfsPerGridFace> scvfsPerGridFace_;
    std::vector<std::vector<GridIndex>> faceIndicesOfScv_;
};


} // end namespace Dumux

#endif
