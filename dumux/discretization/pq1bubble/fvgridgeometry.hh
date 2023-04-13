// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite volume geometry vector for the pq1bubble method
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_GRID_GEOMETRY_HH

#include <utility>
#include <unordered_map>

#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/discretization/method.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/geometry/center.hh>
#include <dumux/geometry/volume.hh>

#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/pq1bubble/geometryhelper.hh>
#include <dumux/discretization/pq1bubble/fvelementgeometry.hh>
#include <dumux/discretization/pq1bubble/subcontrolvolume.hh>
#include <dumux/discretization/pq1bubble/subcontrolvolumeface.hh>
#include <dumux/discretization/pq1bubble/pq1bubblefecache.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

namespace Detail {
template<class GV, class T>
using PQ1BubbleGeometryHelper_t = Dune::Std::detected_or_t<
    Dumux::PQ1BubbleGeometryHelper<GV, typename T::SubControlVolume, typename T::SubControlVolumeFace>,
    SpecifiesGeometryHelper,
    T
>;
} // end namespace Detail

template <class GridView>
struct PQ1BubbleMapperTraits :public DefaultMapperTraits<GridView>
{
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    /**
     * \brief layout for elements and vertices
     *
     */
    static Dune::MCMGLayout layout()
    {
        return [](Dune::GeometryType gt, int dimgrid) {
            return (gt.dim() == dimgrid) || (gt.dim() == 0);
        };
    }
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief The default traits for the pq1bubble finite volume grid geometry
 *        Defines the scv and scvf types and the mapper types
 * \tparam the grid view type
 */
template<class GridView, class MapperTraits = PQ1BubbleMapperTraits<GridView>>
struct PQ1BubbleDefaultGridGeometryTraits
: public MapperTraits
{
    using SubControlVolume = PQ1BubbleSubControlVolume<GridView>;
    using SubControlVolumeFace = PQ1BubbleSubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = PQ1BubbleFVElementGeometry<GridGeometry, enableCache>;
};

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite volume geometry vector for pq1bubble schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class Scalar,
         class GV,
         bool enableCaching = true,
         class Traits = PQ1BubbleDefaultGridGeometryTraits<GV>>
class PQ1BubbleFVGridGeometry
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = PQ1BubbleFVGridGeometry<Scalar, GV, enableCaching, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    static_assert(dim > 1, "Only implemented for dim > 1");

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::PQ1Bubble;
    static constexpr DiscretizationMethod discMethod{};

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::DofMapper;
    //! export the finite element cache type
    using FeCache = Dumux::PQ1BubbleFECache<CoordScalar, Scalar, dim>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor
    PQ1BubbleFVGridGeometry(const GridView gridView)
    : ParentType(gridView)
    , dofMapper_(gridView, Traits::layout())
    , cache_(*this)
    {
        update_();
    }

    //! The dofMapper
    const DofMapper& dofMapper() const
    { return dofMapper_; }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return numScv_; }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->dofMapper().size(); }

    //! update all geometries (call this after grid adaption)
    void update(const GridView& gridView)
    {
        ParentType::update(gridView);
        update_();
    }

    //! update all geometries (call this after grid adaption)
    void update(GridView&& gridView)
    {
        ParentType::update(std::move(gridView));
        update_();
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a vertex / d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicVertexMap_.count(dofIdx); }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicVertexMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries
    const std::unordered_map<GridIndexType, GridIndexType>& periodicVertexMap() const
    { return periodicVertexMap_; }

    //! local view of this object (constructed with the internal cache)
    friend inline LocalView localView(const PQ1BubbleFVGridGeometry& gg)
    { return { gg.cache_ }; }

private:

    class PQ1BubbleGridGeometryCache
    {
        friend class PQ1BubbleFVGridGeometry;
    public:
        //! export the geometry helper type
        using GeometryHelper = Detail::PQ1BubbleGeometryHelper_t<GV, Traits>;

        explicit PQ1BubbleGridGeometryCache(const PQ1BubbleFVGridGeometry& gg)
        : gridGeometry_(&gg)
        {}

        const PQ1BubbleFVGridGeometry& gridGeometry() const
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

        //! Local mappings necessary to construct geometries of scvfs
        const std::vector<std::array<LocalIndexType, 2>>& scvfBoundaryGeometryKeys(GridIndexType eIdx) const
        { return scvfBoundaryGeometryKeys_.at(eIdx); }

    private:
        void clear_()
        {
            scvs_.clear();
            scvfs_.clear();
            hasBoundaryScvf_.clear();
            scvfBoundaryGeometryKeys_.clear();
        }

        std::vector<std::vector<SubControlVolume>> scvs_;
        std::vector<std::vector<SubControlVolumeFace>> scvfs_;
        std::vector<bool> hasBoundaryScvf_;
        std::unordered_map<GridIndexType, std::vector<std::array<LocalIndexType, 2>>> scvfBoundaryGeometryKeys_;

        const PQ1BubbleFVGridGeometry* gridGeometry_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = PQ1BubbleGridGeometryCache;

private:
    using GeometryHelper = typename Cache::GeometryHelper;

    void update_()
    {
        cache_.clear_();
        dofMapper_.update(this->gridView());

        auto numElements = this->gridView().size(0);
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
            auto eIdx = this->elementMapper().index(element);

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            numScv_ += geometryHelper.numScv();
            // construct the sub control volumes
            cache_.scvs_[eIdx].resize(geometryHelper.numScv());
            for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < geometryHelper.numScv(); ++scvLocalIdx)
            {
                auto corners = geometryHelper.getScvCorners(scvLocalIdx);
                cache_.scvs_[eIdx][scvLocalIdx] = SubControlVolume(
                    geometryHelper.scvVolume(scvLocalIdx, corners),
                    geometryHelper.dofPosition(scvLocalIdx),
                    Dumux::center(corners),
                    scvLocalIdx,
                    eIdx,
                    geometryHelper.dofIndex(this->dofMapper(), element, scvLocalIdx),
                    geometryHelper.isOverlappingScv(scvLocalIdx)
                );
            }

            // construct the sub control volume faces
            numScvf_ += geometryHelper.numInteriorScvf();
            cache_.scvfs_[eIdx].resize(geometryHelper.numInteriorScvf());
            LocalIndexType scvfLocalIdx = 0;
            for (; scvfLocalIdx < geometryHelper.numInteriorScvf(); ++scvfLocalIdx)
            {
                const auto scvPair = geometryHelper.getScvPairForScvf(scvfLocalIdx);
                const auto corners = geometryHelper.getScvfCorners(scvfLocalIdx);
                const auto area = Dumux::convexPolytopeVolume<dim-1>(
                    geometryHelper.getInteriorScvfGeometryType(scvfLocalIdx),
                    [&](unsigned int i){ return corners[i]; }
                );

                cache_.scvfs_[eIdx][scvfLocalIdx] = SubControlVolumeFace(
                    Dumux::center(corners),
                    area,
                    geometryHelper.normal(corners, scvPair),
                    std::move(scvPair),
                    scvfLocalIdx,
                    geometryHelper.isOverlappingScvf(scvfLocalIdx)
                );
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
                    cache_.hasBoundaryScvf_[eIdx] = true;

                    const auto localFacetIndex = intersection.indexInInside();
                    const auto numBoundaryScvf = geometryHelper.numBoundaryScvf(localFacetIndex);
                    numScvf_ += numBoundaryScvf;
                    numBoundaryScvf_ += numBoundaryScvf;

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < numBoundaryScvf; ++isScvfLocalIdx)
                    {
                        // find the scvs this scvf is belonging to
                        const auto scvPair = geometryHelper.getScvPairForBoundaryScvf(localFacetIndex, isScvfLocalIdx);
                        const auto corners = geometryHelper.getBoundaryScvfCorners(localFacetIndex, isScvfLocalIdx);
                        const auto area = Dumux::convexPolytopeVolume<dim-1>(
                            Dune::GeometryTypes::cube(dim-1),
                            [&](unsigned int i){ return corners[i]; }
                        );
                        cache_.scvfs_[eIdx].emplace_back(
                            Dumux::center(corners),
                            area,
                            intersection.centerUnitOuterNormal(),
                            std::move(scvPair),
                            scvfLocalIdx,
                            typename SubControlVolumeFace::Traits::BoundaryFlag{ intersection }
                        );

                        // store look-up map to construct boundary scvf geometries
                        cache_.scvfBoundaryGeometryKeys_[eIdx].emplace_back(std::array<LocalIndexType, 2>{{
                            static_cast<LocalIndexType>(localFacetIndex),
                            static_cast<LocalIndexType>(isScvfLocalIdx)
                        }});

                        // increment local counter
                        scvfLocalIdx++;
                    }

                    // TODO also move this to helper class

                    // add all vertices on the intersection to the set of boundary vertices
                    for (int localVIdx = 0; localVIdx < numBoundaryScvf; ++localVIdx)
                    {
                        const auto vIdx = refElement.subEntity(localFacetIndex, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->dofMapper().subIndex(element, vIdx, dim);
                        boundaryDofIndices_[vIdxGlobal] = true;
                    }
                }

                // inform the grid geometry if we have periodic boundaries
                else if (intersection.boundary() && intersection.neighbor())
                {
                    this->setPeriodic();

                    // find the mapped periodic vertex of all vertices on periodic boundaries
                    const auto fIdx = intersection.indexInInside();
                    const auto numFaceVerts = refElement.size(fIdx, 1, dim);
                    const auto eps = 1e-7*(elementGeometry.corner(1) - elementGeometry.corner(0)).two_norm();
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = refElement.subEntity(fIdx, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->dofMapper().subIndex(element, vIdx, dim);
                        const auto vPos = elementGeometry.corner(vIdx);

                        const auto& outside = intersection.outside();
                        const auto outsideGeometry = outside.geometry();
                        for (const auto& isOutside : intersections(this->gridView(), outside))
                        {
                            // only check periodic vertices of the periodic neighbor
                            if (isOutside.boundary() && isOutside.neighbor())
                            {
                                const auto fIdxOutside = isOutside.indexInInside();
                                const auto numFaceVertsOutside = refElement.size(fIdxOutside, 1, dim);
                                for (int localVIdxOutside = 0; localVIdxOutside < numFaceVertsOutside; ++localVIdxOutside)
                                {
                                    const auto vIdxOutside = refElement.subEntity(fIdxOutside, 1, localVIdxOutside, dim);
                                    const auto vPosOutside = outsideGeometry.corner(vIdxOutside);
                                    const auto shift = std::abs((this->bBoxMax()-this->bBoxMin())*intersection.centerUnitOuterNormal());
                                    if (std::abs((vPosOutside-vPos).two_norm() - shift) < eps)
                                        periodicVertexMap_[vIdxGlobal] = this->dofMapper().subIndex(outside, vIdxOutside, dim);
                                }
                            }
                        }
                    }
                }
            }
        }

        // error check: periodic boundaries currently don't work for pq1bubble in parallel
        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for pq1bubble method for parallel simulations!");
    }

    DofMapper dofMapper_;

    const FeCache feCache_;

    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on the boundary
    std::vector<bool> boundaryDofIndices_;

    // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicVertexMap_;

    Cache cache_;
};

} // end namespace Dumux

#endif
