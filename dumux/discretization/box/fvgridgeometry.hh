// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GRID_FVGEOMETRY_HH
#define DUMUX_DISCRETIZATION_BOX_GRID_FVGEOMETRY_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/basefvgridgeometry.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/box/fvelementgeometry.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableFVGridGeometryCache>
class BoxFVGridGeometry
{};

// specialization in case the FVElementGeometries are stored
template<class TypeTag>
class BoxFVGridGeometry<TypeTag, true> : public BaseFVGridGeometry<TypeTag>
{
    using ParentType = BaseFVGridGeometry<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using LocalIndexType = typename SubControlVolumeFace::Traits::LocalIndexType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using Element = typename GridView::template Codim<0>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

    using GeometryHelper = BoxGeometryHelper<GridView, dim, SubControlVolume, SubControlVolumeFace>;

public:
    //! export discretization method
    static constexpr DiscretizationMethods discretizationMethod = DiscretizationMethods::Box;

    //! Constructor
    BoxFVGridGeometry(const GridView gridView)
    : ParentType(gridView) {}

    //! the vertex mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const VertexMapper& dofMapper() const
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
    void update()
    {
        ParentType::update();

        scvs_.clear();
        scvfs_.clear();

        auto numElements = this->gridView().size(0);
        scvs_.resize(numElements);
        scvfs_.resize(numElements);

        boundaryDofIndices_.resize(numDofs(), false);

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
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
            const auto referenceElement = ReferenceElements::general(elementGeometry.type());
#else
            const auto& referenceElement = ReferenceElements::general(elementGeometry.type());
#endif

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            // construct the sub control volumes
            scvs_[eIdx].resize(elementGeometry.corners());
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
                std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(referenceElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                             static_cast<LocalIndexType>(referenceElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

                scvfs_[eIdx][scvfLocalIdx] = SubControlVolumeFace(geometryHelper,
                                                                  element,
                                                                  elementGeometry,
                                                                  scvfLocalIdx,
                                                                  std::move(localScvIndices),
                                                                  false);
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary())
                {
                    const auto isGeometry = intersection.geometry();
                    // count
                    numScvf_ += isGeometry.corners();
                    numBoundaryScvf_ += isGeometry.corners();

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < isGeometry.corners(); ++isScvfLocalIdx)
                    {
                        // find the scvs this scvf is belonging to
                        const LocalIndexType insideScvIdx = static_cast<LocalIndexType>(referenceElement.subEntity(intersection.indexInInside(), 1, isScvfLocalIdx, dim));
                        std::vector<LocalIndexType> localScvIndices = {insideScvIdx, insideScvIdx};

                        scvfs_[eIdx].emplace_back(geometryHelper,
                                                  intersection,
                                                  isGeometry,
                                                  isScvfLocalIdx,
                                                  scvfLocalIdx,
                                                  std::move(localScvIndices),
                                                  true);

                        // increment local counter
                        scvfLocalIdx++;
                    }

                    // add all vertices on the intersection to the set of
                    // boundary vertices
                    const auto fIdx = intersection.indexInInside();
                    const auto numFaceVerts = referenceElement.size(fIdx, 1, dim);
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = referenceElement.subEntity(fIdx, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdx, dim);
                        boundaryDofIndices_[vIdxGlobal] = true;
                    }
                }
            }
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! Get the local scvs for an element
    const std::vector<SubControlVolume>& scvs(IndexType eIdx) const
    { return scvs_[eIdx]; }

    //! Get the local scvfs for an element
    const std::vector<SubControlVolumeFace>& scvfs(IndexType eIdx) const
    { return scvfs_[eIdx]; }

    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(unsigned int dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

private:

    const FeCache feCache_;

    std::vector<std::vector<SubControlVolume>> scvs_;
    std::vector<std::vector<SubControlVolumeFace>> scvfs_;
    // TODO do we need those?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on the boudary
    std::vector<bool> boundaryDofIndices_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class BoxFVGridGeometry<TypeTag, false> : public BaseFVGridGeometry<TypeTag>
{
    using ParentType = BaseFVGridGeometry<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    //! export discretization method
    static constexpr DiscretizationMethods discretizationMethod = DiscretizationMethods::Box;

    //! Constructor
    BoxFVGridGeometry(const GridView gridView)
    : ParentType(gridView) {}

    //! the vertex mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const VertexMapper& dofMapper() const
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
    void update()
    {
        ParentType::update();

        boundaryDofIndices_.resize(numDofs(), false);

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
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
            const auto referenceElement = ReferenceElements::general(elementGeometry.type());
#else
            const auto& referenceElement = ReferenceElements::general(elementGeometry.type());
#endif

            // store the sub control volume face indices on the domain boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                if (intersection.boundary())
                {
                    const auto isGeometry = intersection.geometry();
                    numScvf_ += isGeometry.corners();
                    numBoundaryScvf_ += isGeometry.corners();

                    // add all vertices on the intersection to the set of
                    // boundary vertices
                    const auto fIdx = intersection.indexInInside();
                    const auto numFaceVerts = referenceElement.size(fIdx, 1, dim);
                    for (int localVIdx = 0; localVIdx < numFaceVerts; ++localVIdx)
                    {
                        const auto vIdx = referenceElement.subEntity(fIdx, 1, localVIdx, dim);
                        const auto vIdxGlobal = this->vertexMapper().subIndex(element, vIdx, dim);
                        boundaryDofIndices_[vIdxGlobal] = true;
                    }
                }
            }
        }
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

    //! If a vertex / d.o.f. is on the boundary
    bool dofOnBoundary(unsigned int dofIdx) const
    { return boundaryDofIndices_[dofIdx]; }

private:

    const FeCache feCache_;

    // Information on the global number of geometries
    // TODO do we need those information?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vertices on the boudary
    std::vector<bool> boundaryDofIndices_;
};

} // end namespace Dumux

#endif
