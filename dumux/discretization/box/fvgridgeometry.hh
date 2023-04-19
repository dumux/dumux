// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDiscretization
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GRID_FVGEOMETRY_HH
#define DUMUX_DISCRETIZATION_BOX_GRID_FVGEOMETRY_HH

#include <utility>
#include <unordered_map>
#include <array>
#include <vector>

#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/box/fvelementgeometry.hh>
#include <dumux/discretization/box/subcontrolvolume.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

namespace Detail {
template<class GV, class T>
using BoxGeometryHelper_t = Dune::Std::detected_or_t<
    Dumux::BoxGeometryHelper<GV, GV::dimension, typename T::SubControlVolume, typename T::SubControlVolumeFace>,
    SpecifiesGeometryHelper,
    T
>;
} // end namespace Detail

/*!
 * \ingroup BoxDiscretization
 * \brief The default traits for the box finite volume grid geometry
 *        Defines the scv and scvf types and the mapper types
 * \tparam the grid view type
 */
template<class GridView, class MapperTraits = DefaultMapperTraits<GridView>>
struct BoxDefaultGridGeometryTraits
: public MapperTraits
{
    using SubControlVolume = BoxSubControlVolume<GridView>;
    using SubControlVolumeFace = BoxSubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = BoxFVElementGeometry<GridGeometry, enableCache>;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Base class for the finite volume geometry vector for box schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class Scalar,
         class GridView,
         bool enableGridGeometryCache = false,
         class Traits = BoxDefaultGridGeometryTraits<GridView> >
class BoxFVGridGeometry;

/*!
 * \ingroup BoxDiscretization
 * \brief Base class for the finite volume geometry vector for box schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class Scalar, class GV, class Traits>
class BoxFVGridGeometry<Scalar, GV, true, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = BoxFVGridGeometry<Scalar, GV, true, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::Box;
    static constexpr DiscretizationMethod discMethod{};

    //! export basic grid geometry type for the alternative constructor
    using BasicGridGeometry = BasicGridGeometry_t<GV, Traits>;
    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::VertexMapper;
    //! export the finite element cache type
    using FeCache = Dune::LagrangeLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor with basic grid geometry used to share state with another grid geometry on the same grid view
    BoxFVGridGeometry(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , cache_(*this)
    {
        update_();
    }

    //! Constructor
    BoxFVGridGeometry(const GridView& gridView)
    : BoxFVGridGeometry(std::make_shared<BasicGridGeometry>(gridView))
    {}

    //! the vertex mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const DofMapper& dofMapper() const
    { return this->vertexMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {  return numScv_; }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    //! For compatibility reasons with cc methods
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->vertexMapper().size(); }


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
    friend inline LocalView localView(const BoxFVGridGeometry& gg)
    { return { gg.cache_ }; }

private:

    class BoxGridGeometryCache
    {
        friend class BoxFVGridGeometry;
    public:
        //! export the geometry helper type
        using GeometryHelper = Detail::BoxGeometryHelper_t<GV, Traits>;

        explicit BoxGridGeometryCache(const BoxFVGridGeometry& gg)
        : gridGeometry_(&gg)
        {}

        const BoxFVGridGeometry& gridGeometry() const
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

        //! Returns local mappings for constructing boundary scvf geometries
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

        const BoxFVGridGeometry* gridGeometry_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = BoxGridGeometryCache;

private:
    using GeometryHelper = typename Cache::GeometryHelper;

    void update_()
    {
        cache_.clear_();

        const auto numElements = this->gridView().size(0);
        cache_.scvs_.resize(numElements);
        cache_.scvfs_.resize(numElements);
        cache_.hasBoundaryScvf_.resize(numElements, false);

        boundaryDofIndices_.assign(numDofs(), false);

        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        // Build the SCV and SCV faces
        for (const auto& element : elements(this->gridView()))
        {
            // fill the element map with seeds
            const auto eIdx = this->elementMapper().index(element);

            // count
            numScv_ += element.subEntities(dim);
            numScvf_ += element.subEntities(dim-1);

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            // construct the sub control volumes
            cache_.scvs_[eIdx].resize(elementGeometry.corners());
            for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
            {
                const auto dofIdxGlobal = this->vertexMapper().subIndex(element, scvLocalIdx, dim);

                cache_.scvs_[eIdx][scvLocalIdx] = SubControlVolume(
                    geometryHelper.getScvCorners(scvLocalIdx),
                    scvLocalIdx,
                    eIdx,
                    dofIdxGlobal
                );
            }

            // construct the sub control volume faces
            LocalIndexType scvfLocalIdx = 0;
            cache_.scvfs_[eIdx].resize(element.subEntities(dim-1));
            for (; scvfLocalIdx < element.subEntities(dim-1); ++scvfLocalIdx)
            {
                // find the global and local scv indices this scvf is belonging to
                std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                             static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

                const auto& corners = geometryHelper.getScvfCorners(scvfLocalIdx);
                cache_.scvfs_[eIdx][scvfLocalIdx] = SubControlVolumeFace(
                    corners,
                    geometryHelper.normal(corners, localScvIndices),
                    element,
                    elementGeometry,
                    scvfLocalIdx,
                    std::move(localScvIndices),
                    false
                );
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
                    const auto isGeometry = intersection.geometry();
                    cache_.hasBoundaryScvf_[eIdx] = true;

                    // count
                    numScvf_ += isGeometry.corners();
                    numBoundaryScvf_ += isGeometry.corners();

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < isGeometry.corners(); ++isScvfLocalIdx)
                    {
                        // find the scvs this scvf is belonging to
                        const LocalIndexType insideScvIdx = static_cast<LocalIndexType>(refElement.subEntity(intersection.indexInInside(), 1, isScvfLocalIdx, dim));
                        std::vector<LocalIndexType> localScvIndices = {insideScvIdx, insideScvIdx};

                        cache_.scvfs_[eIdx].emplace_back(
                            geometryHelper.getBoundaryScvfCorners(intersection.indexInInside(), isScvfLocalIdx),
                            intersection.centerUnitOuterNormal(),
                            intersection,
                            isGeometry,
                            isScvfLocalIdx,
                            scvfLocalIdx,
                            std::move(localScvIndices),
                            true
                        );

                        cache_.scvfBoundaryGeometryKeys_[eIdx].emplace_back(std::array<LocalIndexType, 2>{{
                            static_cast<LocalIndexType>(intersection.indexInInside()),
                            static_cast<LocalIndexType>(isScvfLocalIdx)
                        }});

                        // increment local counter
                        scvfLocalIdx++;
                    }

                    // add all vertices on the intersection to the set of
                    // boundary vertices
                    const auto fIdx = intersection.indexInInside();
                    const auto numFaceVerts = refElement.size(fIdx, 1, dim);
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = refElement.subEntity(fIdx, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdx, dim);
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
                        const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdx, dim);
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
                                        periodicVertexMap_[vIdxGlobal] = this->vertexMapper().subIndex(outside, vIdxOutside, dim);
                                }
                            }
                        }
                    }
                }
            }
        }

        // error check: periodic boundaries currently don't work for box in parallel
        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for box method for parallel simulations!");
    }

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

/*!
 * \ingroup BoxDiscretization
 * \brief Base class for the finite volume geometry vector for box schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching disabled we store only some essential index maps to build up local systems on-demand in
 *       the corresponding FVElementGeometry
 */
template<class Scalar, class GV, class Traits>
class BoxFVGridGeometry<Scalar, GV, false, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = BoxFVGridGeometry<Scalar, GV, false, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;

    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::Box;
    static constexpr DiscretizationMethod discMethod{};

    //! export basic grid geometry type for the alternative constructor
    using BasicGridGeometry = BasicGridGeometry_t<GV, Traits>;
    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, false>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::VertexMapper;
    //! export the finite element cache type
    using FeCache = Dune::LagrangeLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor with basic grid geometry used to share state with another grid geometry on the same grid view
    BoxFVGridGeometry(std::shared_ptr<BasicGridGeometry> gg)
    : ParentType(std::move(gg))
    , cache_(*this)
    {
        update_();
    }

    //! Constructor
    BoxFVGridGeometry(const GridView& gridView)
    : BoxFVGridGeometry(std::make_shared<BasicGridGeometry>(gridView))
    {}

    //! the vertex mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const DofMapper& dofMapper() const
    { return this->vertexMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {  return numScv_; }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    //! For compatibility reasons with cc methods
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->vertexMapper().size(); }


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
    friend inline LocalView localView(const BoxFVGridGeometry& gg)
    { return { gg.cache_ }; }

private:

    class BoxGridGeometryCache
    {
        friend class BoxFVGridGeometry;
    public:
        //! export the geometry helper type
        using GeometryHelper = Detail::BoxGeometryHelper_t<GV, Traits>;

        explicit BoxGridGeometryCache(const BoxFVGridGeometry& gg)
        : gridGeometry_(&gg)
        {}

        const BoxFVGridGeometry& gridGeometry() const
        { return *gridGeometry_; }

    private:
        const BoxFVGridGeometry* gridGeometry_;
    };

public:
    //! the cache type (only the caching implementation has this)
    //! this alias should only be used by the local view implementation
    using Cache = BoxGridGeometryCache;

private:

    void update_()
    {
        boundaryDofIndices_.assign(numDofs(), false);

        // save global data on the grid's scvs and scvfs
        // TODO do we need those information?
        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            numScv_ += element.subEntities(dim);
            numScvf_ += element.subEntities(dim-1);

            const auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);

            // store the sub control volume face indices on the domain boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary() && !intersection.neighbor())
                {
                    const auto isGeometry = intersection.geometry();
                    numScvf_ += isGeometry.corners();
                    numBoundaryScvf_ += isGeometry.corners();

                    // add all vertices on the intersection to the set of
                    // boundary vertices
                    const auto fIdx = intersection.indexInInside();
                    const auto numFaceVerts = refElement.size(fIdx, 1, dim);
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = refElement.subEntity(fIdx, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdx, dim);
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
                        const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdx, dim);
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
                                        periodicVertexMap_[vIdxGlobal] = this->vertexMapper().subIndex(outside, vIdxOutside, dim);
                                }
                            }
                        }
                    }
                }
            }
        }

        // error check: periodic boundaries currently don't work for box in parallel
        if (this->isPeriodic() && this->gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Periodic boundaries for box method for parallel simulations!");
    }

    const FeCache feCache_;

    // Information on the global number of geometries
    // TODO do we need those information?
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
