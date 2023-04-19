// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_FV_GRID_GEOMETRY
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_FV_GRID_GEOMETRY

#include <memory>
#include <unordered_map>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/defaultmappertraits.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/geometry/center.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/nonconformingfecache.hh>

#include <dumux/discretization/facecentered/diamond/subcontrolvolume.hh>
#include <dumux/discretization/facecentered/diamond/subcontrolvolumeface.hh>
#include <dumux/discretization/facecentered/diamond/fvelementgeometry.hh>
#include <dumux/discretization/facecentered/diamond/geometryhelper.hh>

namespace Dumux {

namespace Detail {
template<class GV, class T>
using FaceCenteredDiamondGeometryHelper_t = Dune::Std::detected_or_t<
    Dumux::DiamondGeometryHelper<GV, typename T::SubControlVolume, typename T::SubControlVolumeFace>,
    SpecifiesGeometryHelper,
    T
>;
} // end namespace Detail

/*!
 * \ingroup DiamondDiscretization
 * \brief The default traits for the face-centered diamond finite volume grid geometry
 *        Defines the scv and scvf types and the mapper types
 * \tparam GridView the grid view type
 */
template<class GridView>
struct FaceCenteredDiamondDefaultGridGeometryTraits : public DefaultMapperTraits<GridView>
{
    using SubControlVolume = FaceCenteredDiamondSubControlVolume<GridView>;
    using SubControlVolumeFace = FaceCenteredDiamondSubControlVolumeFace<GridView>;
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = FaceCenteredDiamondFVElementGeometry<GridGeometry, enableCache>;
};

/*!
 * \ingroup DiamondDiscretization
 * \brief Grid geometry for the diamond discretization
 */
template<class GV,
         bool enableCaching = true,
         class Traits = FaceCenteredDiamondDefaultGridGeometryTraits<GV>>
class FaceCenteredDiamondFVGridGeometry
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = FaceCenteredDiamondFVGridGeometry<GV, enableCaching, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::SmallLocalIndex;
    using Element = typename GV::template Codim<0>::Entity;

    using Scalar = typename GV::ctype;

    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    static_assert(dim > 1, "Only implemented for dim > 1");

public:
    //! export discretization method
    using DiscretizationMethod = DiscretizationMethods::FCDiamond;
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod{};
    static constexpr bool cachingEnabled = true;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the grid view type
    using GridView = GV;
    //! export the dof mapper type
    using DofMapper = typename Traits::DofMapper;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export the finite element cache type
    using FeCache = NonconformingFECache<Scalar, Scalar, dim>;

    //! Constructor
    FaceCenteredDiamondFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
    : ParentType(gridView)
    , dofMapper_(gridView, Dune::mcmgLayout(Dune::Codim<1>{}))
    , cache_(*this)
    {
        update_();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return numScv_; }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! the total number of dofs
    std::size_t numDofs() const
    { return this->gridView().size(1); }

    //! update all fvElementGeometries (call this after grid adaption)
    void update(const GridView& gridView)
    {
        ParentType::update(gridView);
        update_();
    }

    //! update all fvElementGeometries (call this after grid adaption)
    void update(GridView&& gridView)
    {
        ParentType::update(std::move(gridView));
        update_();
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! If a face / d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! Return a reference to the dof mapper
    const DofMapper& dofMapper() const
    { return dofMapper_; }

    //! If a d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicFaceMap_.count(dofIdx); }

    //! The index of the d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicFaceMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries // TODO rename to periodic dof map in fvassembler
    const std::unordered_map<GridIndexType, GridIndexType>& periodicVertexMap() const
    { return periodicFaceMap_; }

    //! local view of this object (constructed with the internal cache)
    friend inline LocalView localView(const FaceCenteredDiamondFVGridGeometry& gg)
    { return { gg.cache_ }; }

private:

    class FCDiamondGridGeometryCache
    {
        friend class FaceCenteredDiamondFVGridGeometry;
    public:
        //! export the geometry helper type
        using GeometryHelper = Detail::FaceCenteredDiamondGeometryHelper_t<GV, Traits>;

        explicit FCDiamondGridGeometryCache(const FaceCenteredDiamondFVGridGeometry& gg)
        : gridGeometry_(&gg)
        {}

        const FaceCenteredDiamondFVGridGeometry& gridGeometry() const
        { return *gridGeometry_; }

        //! Get the global sub control volume indices of an element
        const std::vector<SubControlVolume>& scvs(GridIndexType eIdx) const
        { return scvs_[eIdx]; }

        //! Get the global sub control volume face indices of an element
        const std::vector<SubControlVolumeFace>& scvfs(GridIndexType eIdx) const
        { return scvfs_[eIdx]; }

        //! Returns whether one of the geometry's scvfs lies on a boundary
        bool hasBoundaryScvf(GridIndexType eIdx) const
        { return hasBoundaryScvf_[eIdx]; }

    private:
        void clear_()
        {
            scvs_.clear();
            scvfs_.clear();
            hasBoundaryScvf_.clear();
        }

        std::vector<std::vector<SubControlVolume>> scvs_;
        std::vector<std::vector<SubControlVolumeFace>> scvfs_;
        std::vector<bool> hasBoundaryScvf_;

        const FaceCenteredDiamondFVGridGeometry* gridGeometry_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = FCDiamondGridGeometryCache;
private:
    using GeometryHelper = typename Cache::GeometryHelper;

    //! update all fvElementGeometries
    void update_()
    {
        // clear containers (necessary after grid refinement)
        cache_.clear_();
        dofMapper_.update(this->gridView());

        // determine size of containers
        const auto numElements = this->gridView().size(0);
        cache_.scvs_.resize(numElements);
        cache_.scvfs_.resize(numElements);
        cache_.hasBoundaryScvf_.resize(numElements, false);

        boundaryDofIndices_.assign(numDofs(), false);

        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;

        // Build the scvs and scv faces
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);

            const auto geometry = element.geometry();
            GeometryHelper geometryHelper(geometry);

            // build the scvs
            cache_.scvs_[eIdx].reserve(geometryHelper.numScv());
            numScv_ += geometryHelper.numScv();
            for (LocalIndexType localScvIdx = 0; localScvIdx < geometryHelper.numScv(); ++localScvIdx)
            {
                const auto dofIndex = dofMapper().subIndex(element, localScvIdx, 1);
                const auto& corners = geometryHelper.getScvCorners(localScvIdx);
                const auto volume = Dumux::convexPolytopeVolume<dim>(
                    SubControlVolume::Traits::geometryType(geometry.type()),
                    [&](unsigned int i){ return corners[i]; }
                );

                cache_.scvs_[eIdx].emplace_back(
                    volume,
                    geometryHelper.facetCenter(localScvIdx),
                    Dumux::center(corners),
                    localScvIdx,
                    eIdx,
                    dofIndex
                );
            }

            // build interior scvfs
            LocalIndexType localScvfIdx = 0;
            cache_.scvfs_[eIdx].reserve(geometryHelper.numInteriorScvf());
            numScvf_ += geometryHelper.numInteriorScvf();
            for (; localScvfIdx < geometryHelper.numInteriorScvf(); ++localScvfIdx)
            {
                const auto& corners = geometryHelper.getScvfCorners(localScvfIdx);
                const auto& scvPair = geometryHelper.getInsideOutsideScvForScvf(localScvfIdx);
                const auto area = Dumux::convexPolytopeVolume<dim-1>(
                    SubControlVolumeFace::Traits::interiorGeometryType(geometry.type()),
                    [&](unsigned int i){ return corners[i]; }
                );

                cache_.scvfs_[eIdx].emplace_back(
                    Dumux::center(corners),
                    area,
                    geometryHelper.normal(corners, scvPair),
                    scvPair,
                    localScvfIdx
                );
            }

            // build boundary scvfs
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (onDomainBoundary_(intersection))
                {
                    // store information that the face dof is on a boundary
                    const LocalIndexType localFacetIndex = intersection.indexInInside();
                    const auto dofIndex = dofMapper().subIndex(element, localFacetIndex, 1);
                    boundaryDofIndices_[dofIndex] = true;

                    // and that the element has a boundary face
                    cache_.hasBoundaryScvf_[eIdx] = true;

                    // add boundary scvf
                    const auto geo = intersection.geometry();
                    cache_.scvfs_[eIdx].emplace_back(
                        geo.center(),
                        geo.volume(),
                        intersection.centerUnitOuterNormal(),
                        std::array<LocalIndexType, 2>{{localFacetIndex, localFacetIndex}},
                        localScvfIdx,
                        typename SubControlVolumeFace::Traits::BoundaryFlag{ intersection }
                    );

                    // increment local and global counters
                    ++localScvfIdx;
                    ++numBoundaryScvf_;
                    ++numScvf_;
                }

                // handle periodic boundaries
                if (onPeriodicBoundary_(intersection))
                {
                    const LocalIndexType localFacetIndex = intersection.indexInInside();
                    const auto dofIndex = dofMapper().subIndex(element, localFacetIndex, 1);

                    this->setPeriodic();

                    const auto& otherElement = intersection.outside();

                    LocalIndexType otherIntersectionLocalIdx = 0;
                    bool periodicFaceFound = false;

                    for (const auto& otherIntersection : intersections(this->gridView(), otherElement))
                    {
                        if (periodicFaceFound)
                            continue;

                        if (Dune::FloatCmp::eq(intersection.centerUnitOuterNormal()*otherIntersection.centerUnitOuterNormal(), -1.0, 1e-7))
                        {
                            const auto periodicDofIdx = dofMapper().subIndex(otherElement, otherIntersectionLocalIdx, 1);
                            periodicFaceMap_[dofIndex] = periodicDofIdx;
                            periodicFaceFound = true;
                        }

                        ++otherIntersectionLocalIdx;
                    }
                }
            }
        }
    }

    bool onDomainBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && intersection.boundary();
    }

    bool onProcessorBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && !intersection.boundary();
    }

    bool onPeriodicBoundary_(const typename GridView::Intersection& intersection) const
    {
        return intersection.boundary() && intersection.neighbor();
    }

    // faces on the boundary
    std::vector<bool> boundaryDofIndices_;

    DofMapper dofMapper_;

    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicFaceMap_;

    const FeCache feCache_;

    Cache cache_;
};

} // end namespace Dumux

#endif
