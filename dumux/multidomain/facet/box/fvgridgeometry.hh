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
 * \ingroup FacetCoupling
 * \brief Base class for the finite volume grid geometry for box models in the
 *        context of models considering coupling of different domains across the
 *        bulk grid facets. This builds up the sub control volumes and sub control
 *        volume faces for each element of the grid partition.
 */
#ifndef DUMUX_FACETCOUPLING_BOX_GRID_FVGEOMETRY_HH
#define DUMUX_FACETCOUPLING_BOX_GRID_FVGEOMETRY_HH

#include <algorithm>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/box/subcontrolvolume.hh>

#include <dumux/multidomain/facet/box/fvelementgeometry.hh>
#include <dumux/multidomain/facet/box/subcontrolvolumeface.hh>
#include <dumux/multidomain/facet/vertexmapper.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief The default traits for the finite volume grid geometry
 *        of the box scheme with coupling occuring across the element facets.
 *        Defines the scv and scvf types and the mapper types.
 * \tparam the grid view type
 */
template<class GridView>
struct BoxFacetCouplingDefaultGridGeometryTraits
{
    // use a specialized version of the box scvf
    using SubControlVolume = BoxSubControlVolume<GridView>;
    using SubControlVolumeFace = BoxFacetCouplingSubControlVolumeFace<GridView>;

    template<class GridGeometry, bool enableCache>
    using LocalView = BoxFacetCouplingFVElementGeometry<GridGeometry, enableCache>;

    // per default we use an mcmg mapper for the elements
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    // the default vertex mapper is the enriched vertex dof mapper
    using VertexMapper = EnrichedVertexDofMapper<GridView>;
    // Mapper type for mapping edges
    using FacetMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
};

/*!
 * \ingroup FacetCoupling
 * \brief Base class for the finite volume geometry vector for box schemes in the context
 *        of coupled models where the coupling occurs across the element facets. This builds
 *        up the sub control volumes and sub control volume faces.
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class Scalar,
         class GridView,
         bool enableGridGeometryCache = false,
         class Traits = BoxFacetCouplingDefaultGridGeometryTraits<GridView> >
class BoxFacetCouplingFVGridGeometry;

/*!
 * \ingroup FacetCoupling
 * \brief Base class for the finite volume geometry vector for box schemes in the context
 *        of coupled models where the coupling occurs across the element facets. This builds
 *        up the sub control volumes and sub control volume faces.
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class Scalar, class GV, class Traits>
class BoxFacetCouplingFVGridGeometry<Scalar, GV, true, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = BoxFacetCouplingFVGridGeometry<Scalar, GV, true, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

    using Element = typename GV::template Codim<0>::Entity;
    using CoordScalar = typename GV::ctype;
    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

    using GeometryHelper = BoxGeometryHelper<GV, dim, typename Traits::SubControlVolume, typename Traits::SubControlVolumeFace>;

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::box;

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
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor
    BoxFacetCouplingFVGridGeometry(const GridView& gridView)
    : ParentType(gridView) {}

    //! the vertex mapper is the dofMapper
    const DofMapper& dofMapper() const
    { return this->vertexMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return numScv_; }

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

    /*!
     * \brief update all fvElementGeometries (do this again after grid adaption)
     * \note This assumes conforming grids!
     *
     * \param facetGridView The grid view of a (dim-1)-dimensional grid conforming
     *                      with the facets of this grid view, indicating on which facets
     *                      nodal dofs should be enriched.
     * \param codimOneGridAdapter Adapter class that allows access to information on the d-
     *                            dimensional grid for entities of the (d-1)-dimensional grid
     * \param verbose Verbosity level for vertex enrichment
     */
    template<class FacetGridView, class CodimOneGridAdapter>
    void update(const FacetGridView& facetGridView,
                const CodimOneGridAdapter& codimOneGridAdapter,
                bool verbose = false)
    {
        // first update the parent (mappers etc)
        ParentType::update();

        // enrich the vertex mapper subject to the provided facet grid
        this->vertexMapper().enrich(facetGridView, codimOneGridAdapter, verbose);

        // resize containers
        const auto numDof = numDofs();
        const auto numElements = this->gridView().size(0);
        scvs_.clear();
        scvfs_.clear();
        scvs_.resize(numElements);
        scvfs_.resize(numElements);
        boundaryDofIndices_.assign(numDof, false);
        interiorBoundaryDofIndices_.assign(numDof, false);

        // Build the SCV and SCV faces
        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            auto eIdx = this->elementMapper().index(element);

            // keep track of number of scvs and scvfs
            numScv_ += element.subEntities(dim);
            numScvf_ += element.subEntities(dim-1);

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            // construct the sub control volumes
            scvs_[eIdx].clear();
            scvs_[eIdx].reserve(elementGeometry.corners());
            for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
                scvs_[eIdx].emplace_back(geometryHelper,
                                         scvLocalIdx,
                                         eIdx,
                                         this->vertexMapper().subIndex(element, scvLocalIdx, dim));

            // construct the sub control volume faces
            LocalIndexType scvfLocalIdx = 0;
            scvfs_[eIdx].clear();
            scvfs_[eIdx].reserve(element.subEntities(dim-1));
            for (; scvfLocalIdx < element.subEntities(dim-1); ++scvfLocalIdx)
            {
                // find the global and local scv indices this scvf is belonging to
                std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                             static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

                // create the sub-control volume face
                scvfs_[eIdx].emplace_back(geometryHelper,
                                          element,
                                          elementGeometry,
                                          scvfLocalIdx,
                                          std::move(localScvIndices));
            }

            // construct the sub control volume faces on the domain/interior boundaries
            // skip handled facets (necessary for e.g. Dune::FoamGrid)
            std::vector<unsigned int> handledFacets;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (std::count(handledFacets.begin(), handledFacets.end(), intersection.indexInInside()))
                    continue;

                handledFacets.push_back(intersection.indexInInside());

                // determine if all corners live on the facet grid
                const auto isGeometry = intersection.geometry();
                const auto numFaceCorners = isGeometry.corners();
                const auto idxInInside = intersection.indexInInside();
                const auto boundary = intersection.boundary();

                std::vector<LocalIndexType> vIndicesLocal(numFaceCorners);
                for (int i = 0; i < numFaceCorners; ++i)
                    vIndicesLocal[i] = static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, i, dim));

                std::vector<LocalIndexType> gridVertexIndices(numFaceCorners);
                for (int i = 0; i < numFaceCorners; ++i)
                    gridVertexIndices[i] = this->vertexMapper().vertexIndex(element, vIndicesLocal[i], dim);

                // if the vertices compose a facet element, the intersection is on facet grid
                const bool isOnFacet = codimOneGridAdapter.composeFacetElement(gridVertexIndices);

                // make sure there are no periodic boundaries
                if (boundary && intersection.neighbor())
                    DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the box facet coupling scheme");

                // if it is not, but it is on the boundary -> boundary scvf
                if (isOnFacet || boundary)
                {
                    // keep track of number of faces
                    numScvf_ += numFaceCorners;
                    numBoundaryScvf_ += int(boundary)*numFaceCorners;

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < numFaceCorners; ++isScvfLocalIdx)
                    {
                        // find the inside scv this scvf is belonging to (localIdx = element local vertex index)
                        std::vector<LocalIndexType> localScvIndices = {vIndicesLocal[isScvfLocalIdx], vIndicesLocal[isScvfLocalIdx]};

                        // create the sub-control volume face
                        scvfs_[eIdx].emplace_back(geometryHelper,
                                                  intersection,
                                                  isGeometry,
                                                  isScvfLocalIdx,
                                                  scvfLocalIdx,
                                                  std::move(localScvIndices),
                                                  boundary,
                                                  isOnFacet);

                        // Mark vertices to be on domain and/or interior boundary
                        const auto dofIndex = this->vertexMapper().subIndex(element, vIndicesLocal[isScvfLocalIdx], dim);
                        if (boundary) boundaryDofIndices_[ dofIndex ] = boundary;
                        if (isOnFacet) interiorBoundaryDofIndices_[ dofIndex ] = isOnFacet;

                        // increment local counter
                        scvfLocalIdx++;
                    }
                }
            }
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! Get the local scvs for an element
    const std::vector<SubControlVolume>& scvs(GridIndexType eIdx) const
    { return scvs_[eIdx]; }

    //! Get the local scvfs for an element
    const std::vector<SubControlVolumeFace>& scvfs(GridIndexType eIdx) const
    { return scvfs_[eIdx]; }

    //! If a d.o.f. is on the boundary
    bool dofOnBoundary(GridIndexType dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a d.o.f. is on an interior boundary
    bool dofOnInteriorBoundary(GridIndexType dofIdx) const
    { return interiorBoundaryDofIndices_[dofIdx]; }

    //! Periodic boundaries are not supported for the box facet coupling scheme
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return false; }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the box facet coupling scheme"); }

    //! Returns the map between dofs across periodic boundaries
    std::unordered_map<GridIndexType, GridIndexType> periodicVertexMap() const
    { return std::unordered_map<GridIndexType, GridIndexType>(); }

private:
    const FeCache feCache_;

    std::vector<std::vector<SubControlVolume>> scvs_;
    std::vector<std::vector<SubControlVolumeFace>> scvfs_;

    // TODO do we need those?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on domain/interior boundaries
    std::vector<bool> boundaryDofIndices_;
    std::vector<bool> interiorBoundaryDofIndices_;
};

/*!
 * \ingroup FacetCoupling
 * \brief Base class for the finite volume geometry vector for box schemes
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching disabled we store only some essential index maps to build up local systems on-demand in
 *       the corresponding FVElementGeometry
 */
template<class Scalar, class GV, class Traits>
class BoxFacetCouplingFVGridGeometry<Scalar, GV, false, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = BoxFacetCouplingFVGridGeometry<Scalar, GV, false, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

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
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export dof mapper type
    using DofMapper = typename Traits::VertexMapper;
    //! export the finite element cache type
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    //! export the grid view type
    using GridView = GV;

    //! Constructor
    BoxFacetCouplingFVGridGeometry(const GridView gridView)
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
    { return this->vertexMapper().size(); }

    /*!
     * \brief update all fvElementGeometries (do this again after grid adaption)
     * \note This assumes conforming grids!
     *
     * \param facetGridView The grid view of a (dim-1)-dimensional grid conforming
     *                      with the facets of this grid view, indicating on which facets
     *                      nodal dofs should be enriched.
     * \param codimOneGridAdapter Adapter class that allows access to information on the d-
     *                            dimensional grid for entities of the (d-1)-dimensional grid
     * \param verbose Verbosity level
     */
    template<class FacetGridView, class CodimOneGridAdapter>
    void update(const FacetGridView& facetGridView,
                const CodimOneGridAdapter& codimOneGridAdapter,
                bool verbose = false)
    {
        // first update the parent (mappers etc)
        ParentType::update();

        // enrich the vertex mapper subject to the provided facet grid
        this->vertexMapper().enrich(facetGridView, codimOneGridAdapter, verbose);

        // save global data on the grid's scvs and scvfs
        // TODO do we need those information?
        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;

        const auto numDof = numDofs();
        boundaryDofIndices_.assign(numDof, false);
        interiorBoundaryDofIndices_.assign(numDof, false);
        facetIsOnInteriorBoundary_.assign(this->gridView().size(1), false);
        for (const auto& element : elements(this->gridView()))
        {
            numScv_ += element.subEntities(dim);
            numScvf_ += element.subEntities(dim-1);

            const auto elementGeometry = element.geometry();
            const auto refElement = referenceElement(elementGeometry);

            // store the sub control volume face indices on the domain/interior boundary
            // skip handled facets (necessary for e.g. Dune::FoamGrid)
            std::vector<unsigned int> handledFacets;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (std::count(handledFacets.begin(), handledFacets.end(), intersection.indexInInside()))
                    continue;

                handledFacets.push_back(intersection.indexInInside());

                // determine if all corners live on the facet grid
                const auto isGeometry = intersection.geometry();
                const auto numFaceCorners = isGeometry.corners();
                const auto idxInInside = intersection.indexInInside();
                const auto boundary = intersection.boundary();

                std::vector<LocalIndexType> vIndicesLocal(numFaceCorners);
                for (int i = 0; i < numFaceCorners; ++i)
                    vIndicesLocal[i] = static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, i, dim));

                std::vector<GridIndexType> gridVertexIndices(numFaceCorners);
                for (int i = 0; i < numFaceCorners; ++i)
                    gridVertexIndices[i] = this->vertexMapper().vertexIndex(element, vIndicesLocal[i], dim);

                // if all vertices are living on the facet grid, this is an interiour boundary
                const bool isOnFacet = codimOneGridAdapter.composeFacetElement(gridVertexIndices);

                // make sure there are no periodic boundaries
                if (boundary && intersection.neighbor())
                    DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the box facet coupling scheme");

                if (isOnFacet || boundary)
                {
                    numScvf_ += numFaceCorners;
                    numBoundaryScvf_ += int(boundary)*numFaceCorners;

                    // Mark vertices to be on domain and/or interior boundary
                    for (int i = 0; i < numFaceCorners; ++i)
                    {
                        const auto dofIndex = this->vertexMapper().subIndex(element, vIndicesLocal[i], dim);
                        if (boundary) boundaryDofIndices_[ dofIndex ] = true;
                        if (isOnFacet)
                        {
                            interiorBoundaryDofIndices_[ dofIndex ] = true;
                            facetIsOnInteriorBoundary_[ facetMapper_.subIndex(element, idxInInside, 1) ] = true;
                        }
                    }
                }
            }
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! If a d.o.f. is on the boundary
    bool dofOnBoundary(unsigned int dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

    //! If a d.o.f. is on an interior boundary
    bool dofOnInteriorBoundary(unsigned int dofIdx) const
    { return interiorBoundaryDofIndices_[dofIdx]; }

    //! returns true if an intersection is on an interior boundary
    bool isOnInteriorBoundary(const Element& element, const Intersection& intersection) const
    { return facetIsOnInteriorBoundary_[ facetMapper_.subIndex(element, intersection.indexInInside(), 1) ]; }

    //! Periodic boundaries are not supported for the box facet coupling scheme
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return false; }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Periodic boundaries are not supported by the facet coupling scheme"); }

    //! Returns the map between dofs across periodic boundaries
    std::unordered_map<GridIndexType, GridIndexType> periodicVertexMap() const
    { return std::unordered_map<GridIndexType, GridIndexType>(); }

private:
    const FeCache feCache_;

    // Information on the global number of geometries
    // TODO do we need those information?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on domain/interior boundaries
    std::vector<bool> boundaryDofIndices_;
    std::vector<bool> interiorBoundaryDofIndices_;

    // facet mapper and markers which facets lie on interior boundaries
    typename Traits::FacetMapper facetMapper_;
    std::vector<bool> facetIsOnInteriorBoundary_;
};

} // end namespace Dumux

#endif
