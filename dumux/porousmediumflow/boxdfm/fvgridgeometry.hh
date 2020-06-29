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
 * \ingroup BoxDFMModel
 * \brief Base class for the finite volume geometry vector for box schemes that consider
 *        extra connectivity between grid vertices on marked codim one entities.
 *
 * On these, an additional scvf is created accounting for the additional exchange fluxes
 * between these degrees of freedom.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_GRID_FVGEOMETRY_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_GRID_FVGEOMETRY_HH

#include <unordered_map>

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/discretization/method.hh>
#include <dumux/common/defaultmappertraits.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/extrusion.hh>

#include "fvelementgeometry.hh"
#include "geometryhelper.hh"
#include "subcontrolvolume.hh"
#include "subcontrolvolumeface.hh"

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief The default traits for the box finite volume grid geometry
 *
 * Defines the scv and scvf types and the mapper types.
 *
 * \tparam the grid view type
 */
template<class GridView, class MapperTraits = DefaultMapperTraits<GridView>>
struct BoxDfmDefaultGridGeometryTraits
: public MapperTraits
{
    using SubControlVolume = BoxDfmSubControlVolume<GridView>;
    using SubControlVolumeFace = BoxDfmSubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = BoxDfmFVElementGeometry<GridGeometry, enableCache>;

    // Mapper type for mapping edges
    using FacetMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
};

/*!
 * \ingroup BoxDFMModel
 * \brief Base class for the finite volume geometry vector for box schemes
 *
 * This builds up the sub control volumes and sub control volume faces
 *
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class Scalar,
         class GridView,
         bool enableGridGeometryCache = false,
         class Traits = BoxDfmDefaultGridGeometryTraits<GridView> >
class BoxDfmFVGridGeometry;

/*!
 * \ingroup BoxDFMModel
 * \brief Base class for the finite volume geometry vector for box schemes that consider
 *        extra connectivity between grid vertices on marked codim one entities.
 *
 * On these, an additional scvf is created accounting for the additional exchange fluxes
 * between these degrees of freedom.
 *
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class Scalar, class GV, class Traits>
class BoxDfmFVGridGeometry<Scalar, GV, true, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = BoxDfmFVGridGeometry<Scalar, GV, true, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename GV::IndexSet::IndexType;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;
    static_assert(dim == 2 || dim == 3, "The box-dfm GridGeometry is only implemented in 2 or 3 dimensions.");

    using GeometryHelper = BoxDfmGeometryHelper<GV, dim,
                                                typename Traits::SubControlVolume,
                                                typename Traits::SubControlVolumeFace>;

public:
    //! Export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::box;

    //! Export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! Export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! Export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! Export the extrusion type
    using Extrusion = Extrusion_t<Traits>;
    //! Export dof mapper type
    using DofMapper = typename Traits::VertexMapper;
    //! Export the finite element cache type
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! Export the grid view type
    using GridView = GV;

    //! Constructor
    BoxDfmFVGridGeometry(const GridView gridView)
    : ParentType(gridView) {}

    //! The vertex mapper is the dofMapper
    //! This is convenience to have better chance to have the same main files for box/tpfa/mpfa...
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
    { return this->gridView().size(dim); }

    //! Update all fvElementGeometries (do this again after grid adaption)
    template< class FractureGridAdapter >
    void update(const FractureGridAdapter& fractureGridAdapter)
    {
        ParentType::update();

        scvs_.clear();
        scvfs_.clear();

        auto numElements = this->gridView().size(0);
        scvs_.resize(numElements);
        scvfs_.resize(numElements);

        boundaryDofIndices_.assign(numDofs(), false);
        fractureDofIndices_.assign(this->gridView.size(dim), false);

        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        // Build the SCV and SCV faces
        for (const auto& element : elements(this->gridView()))
        {
            // fill the element map with seeds
            auto eIdx = this->elementMapper().index(element);

            // count
            numScv_ += element.subEntities(dim);
            numScvf_ += element.subEntities(dim-1);

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            // construct the sub control volumes
            scvs_[eIdx].resize(elementGeometry.corners());
            using LocalIndexType = typename SubControlVolumeFace::Traits::LocalIndexType;
            for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
            {
                const auto dofIdxGlobal = this->vertexMapper().subIndex(element, scvLocalIdx, dim);

                scvs_[eIdx][scvLocalIdx] = SubControlVolume(geometryHelper,
                                                            scvLocalIdx,
                                                            eIdx,
                                                            dofIdxGlobal);
            }

            // construct the sub control volume faces
            LocalIndexType scvfLocalIdx = 0;
            scvfs_[eIdx].resize(element.subEntities(dim-1));
            for (; scvfLocalIdx < element.subEntities(dim-1); ++scvfLocalIdx)
            {
                // find the global and local scv indices this scvf is belonging to
                std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                             static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

                scvfs_[eIdx][scvfLocalIdx] = SubControlVolumeFace(geometryHelper,
                                                                  element,
                                                                  elementGeometry,
                                                                  scvfLocalIdx,
                                                                  std::move(localScvIndices));
            }

            // construct the ...
            // ... sub-control volume faces on the domain boundary
            // ... sub-control volumes on fracture facets
            // ... sub-control volume faces on fracture facets
            // NOTE We do not construct fracture scvfs on boundaries here!
            //      That means specifying Neumann fluxes on only the fractures is not possible
            //      However, it is difficult to interpret this here as a fracture ending on the boundary
            //      could also be connected to a facet which is both boundary and fracture at the same time!
            //      In that case, the fracture boundary scvf wouldn't make sense. In order to do it properly
            //      we would have to find only those fractures that are at the boundary and aren't connected
            //      to a fracture which is a boundary.
            LocalIndexType scvLocalIdx = element.subEntities(dim);
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                // first, obtain all vertex indices on this intersection
                const auto& isGeometry = intersection.geometry();
                const auto numCorners = isGeometry.corners();
                const auto idxInInside = intersection.indexInInside();

                std::vector<GridIndexType> isVertexIndices(numCorners);
                for (unsigned int vIdxLocal = 0; vIdxLocal < numCorners; ++vIdxLocal)
                    isVertexIndices[vIdxLocal] = this->vertexMapper().subIndex(element,
                                                                               refElement.subEntity(idxInInside, 1, vIdxLocal, dim),
                                                                               dim);
                // maybe add boundary scvf
                if (intersection.boundary() && !intersection.neighbor())
                {
                    numScvf_ += isGeometry.corners();
                    numBoundaryScvf_ += isGeometry.corners();

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < numCorners; ++isScvfLocalIdx)
                    {
                        // find the scvs this scvf is belonging to
                        const LocalIndexType insideScvIdx = static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, isScvfLocalIdx, dim));
                        std::vector<LocalIndexType> localScvIndices = {insideScvIdx, insideScvIdx};
                        scvfs_[eIdx].emplace_back(geometryHelper,
                                                  intersection,
                                                  isGeometry,
                                                  isScvfLocalIdx,
                                                  scvfLocalIdx++,
                                                  std::move(localScvIndices));
                    }

                    // add all vertices on the intersection to the set of
                    // boundary vertices
                    const auto numFaceVerts = refElement.size(idxInInside, 1, dim);
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = refElement.subEntity(idxInInside, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdx, dim);
                        boundaryDofIndices_[vIdxGlobal] = true;
                    }
                }
                // make sure we have no periodic boundaries
                else if (intersection.boundary() && intersection.neighbor())
                    DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the box-dfm scheme");

                // maybe add fracture scvs & scvfs
                if (fractureGridAdapter.composeFacetElement(isVertexIndices))
                {
                    for (auto vIdx : isVertexIndices)
                        fractureDofIndices_[vIdx] = true;

                    // add fracture scv for each vertex of intersection
                    numScv_ += numCorners;
                    const auto curNumScvs = scvs_[eIdx].size();
                    scvs_[eIdx].reserve(curNumScvs+numCorners);
                    for (unsigned int vIdxLocal = 0; vIdxLocal < numCorners; ++vIdxLocal)
                        scvs_[eIdx].emplace_back(geometryHelper,
                                                 intersection,
                                                 isGeometry,
                                                 vIdxLocal,
                                                 static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, vIdxLocal, dim)),
                                                 scvLocalIdx++,
                                                 idxInInside,
                                                 eIdx,
                                                 isVertexIndices[vIdxLocal]);

                    // add fracture scvf for each edge of the intersection in 3d
                    if (dim == 3)
                    {
                        const auto& faceRefElement = referenceElement(isGeometry);
                        for (unsigned int edgeIdx = 0; edgeIdx < faceRefElement.size(1); ++edgeIdx)
                        {
                            // inside/outside scv indices in face local node numbering
                            std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(faceRefElement.subEntity(edgeIdx, 1, 0, dim-1)),
                                                                         static_cast<LocalIndexType>(faceRefElement.subEntity(edgeIdx, 1, 1, dim-1))});

                            // add offset to get the right scv indices
                            std::for_each( localScvIndices.begin(),
                                           localScvIndices.end(),
                                           [curNumScvs] (auto& elemLocalIdx) { elemLocalIdx += curNumScvs; } );

                            // add scvf
                            numScvf_++;
                            scvfs_[eIdx].emplace_back(geometryHelper,
                                                      intersection,
                                                      isGeometry,
                                                      edgeIdx,
                                                      scvfLocalIdx++,
                                                      std::move(localScvIndices),
                                                      intersection.boundary());
                        }
                    }

                    // dim == 2, intersection is an edge, make 1 scvf
                    else
                    {
                        // inside/outside scv indices in face local node numbering
                        std::vector<LocalIndexType> localScvIndices({0, 1});

                        // add offset such that the fracture scvs above are addressed
                        std::for_each( localScvIndices.begin(),
                                       localScvIndices.end(),
                                       [curNumScvs] (auto& elemLocalIdx) { elemLocalIdx += curNumScvs; } );

                        // add scvf
                        numScvf_++;
                        scvfs_[eIdx].emplace_back(geometryHelper,
                                                  intersection,
                                                  isGeometry,
                                                  /*idxOnIntersection*/0,
                                                  scvfLocalIdx++,
                                                  std::move(localScvIndices),
                                                  intersection.boundary());
                    }
                }
            }
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const { return feCache_; }
    //! Get the local scvs for an element
    const std::vector<SubControlVolume>& scvs(GridIndexType eIdx) const { return scvs_[eIdx]; }
    //! Get the local scvfs for an element
    const std::vector<SubControlVolumeFace>& scvfs(GridIndexType eIdx) const { return scvfs_[eIdx]; }
    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(unsigned int dofIdx) const { return boundaryDofIndices_[dofIdx]; }
    //! If a vertex / d.o.f. is on a fracture
    bool dofOnFracture(unsigned int dofIdx) const { return fractureDofIndices_[dofIdx]; }
    //! Periodic boundaries are not supported for the box-dfm scheme
    bool dofOnPeriodicBoundary(std::size_t dofIdx) const { return false; }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    std::size_t periodicallyMappedDof(std::size_t dofIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the box-dfm scheme"); }

    //! Returns the map between dofs across periodic boundaries
    std::unordered_map<std::size_t, std::size_t> periodicVertexMap() const
    { return std::unordered_map<std::size_t, std::size_t>(); }

private:
    const FeCache feCache_;

    std::vector<std::vector<SubControlVolume>> scvs_;
    std::vector<std::vector<SubControlVolumeFace>> scvfs_;

    // TODO do we need those?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on the boundary & fracture facets
    std::vector<bool> boundaryDofIndices_;
    std::vector<bool> fractureDofIndices_;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Base class for the finite volume geometry vector for box schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching disabled we store only some essential index maps to build up local systems on-demand in
 *       the corresponding FVElementGeometry
 */
template<class Scalar, class GV, class Traits>
class BoxDfmFVGridGeometry<Scalar, GV, false, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = BoxDfmFVGridGeometry<Scalar, GV, false, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename GV::IndexSet::IndexType;

    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    using Element = typename GV::template Codim<0>::Entity;
    using Intersection = typename GV::Intersection;
    using CoordScalar = typename GV::ctype;

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::box;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, false>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! Export the extrusion type
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::VertexMapper;
    //! export the finite element cache type
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor
    BoxDfmFVGridGeometry(const GridView gridView)
    : ParentType(gridView)
    , facetMapper_(gridView, Dune::mcmgLayout(Dune::template Codim<1>()))
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
    { return this->gridView().size(dim); }

    //! update all fvElementGeometries (do this again after grid adaption)
    template< class FractureGridAdapter >
    void update(const FractureGridAdapter& fractureGridAdapter)
    {
        ParentType::update();

        boundaryDofIndices_.assign(numDofs(), false);
        fractureDofIndices_.assign(numDofs(), false);
        facetOnFracture_.assign(this->gridView().size(1), false);

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
                // first, obtain all vertex indices on this intersection
                const auto& isGeometry = intersection.geometry();
                const auto numCorners = isGeometry.corners();
                const auto idxInInside = intersection.indexInInside();

                std::vector<GridIndexType> isVertexIndices(numCorners);
                for (unsigned int vIdxLocal = 0; vIdxLocal < numCorners; ++vIdxLocal)
                    isVertexIndices[vIdxLocal] = this->vertexMapper().subIndex(element,
                                                                               refElement.subEntity(idxInInside, 1, vIdxLocal, dim),
                                                                               dim);

                if (intersection.boundary() && !intersection.neighbor())
                {
                    numScvf_ += numCorners;
                    numBoundaryScvf_ += numCorners;

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

                // make sure we have no periodic boundaries
                else if (intersection.boundary() && intersection.neighbor())
                    DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the box-dfm scheme");

                // maybe add fracture scvs & scvfs
                if (fractureGridAdapter.composeFacetElement(isVertexIndices))
                {
                    facetOnFracture_[facetMapper_.subIndex(element, idxInInside, 1)] = true;
                    for (auto vIdx : isVertexIndices)
                        fractureDofIndices_[vIdx] = true;

                    const auto isGeometry = intersection.geometry();
                    numScv_ += isGeometry.corners();
                    numScvf_ += dim == 3 ? referenceElement(isGeometry).size(1) : 1;
                }
            }
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const { return feCache_; }
    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(unsigned int dofIdx) const { return boundaryDofIndices_[dofIdx]; }
    //! If a vertex / d.o.f. is on a fracture
    bool dofOnFracture(unsigned int dofIdx) const { return fractureDofIndices_[dofIdx]; }
    //! Periodic boundaries are not supported for the box-dfm scheme
    bool dofOnPeriodicBoundary(std::size_t dofIdx) const { return false; }

    //! Returns true if an intersection coincides with a fracture element
    bool isOnFracture(const Element& element, const Intersection& intersection) const
    { return facetOnFracture_[facetMapper_.subIndex(element, intersection.indexInInside(), 1)]; }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    std::size_t periodicallyMappedDof(std::size_t dofIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the box-dfm scheme"); }

    //! Returns the map between dofs across periodic boundaries
    std::unordered_map<std::size_t, std::size_t> periodicVertexMap() const
    { return std::unordered_map<std::size_t, std::size_t>(); }

private:

    const FeCache feCache_;

    // Information on the global number of geometries
    // TODO do we need those information?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on the boundary & fracture facets
    std::vector<bool> boundaryDofIndices_;
    std::vector<bool> fractureDofIndices_;

    // facet mapper and markers which facets lie on interior boundaries
    typename Traits::FacetMapper facetMapper_;
    std::vector<bool> facetOnFracture_;
};

} // end namespace Dumux

#endif
