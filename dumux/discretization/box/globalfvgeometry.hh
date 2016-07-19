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
#ifndef DUMUX_DISCRETIZATION_BOX_GLOBAL_FVGEOMETRY_HH
#define DUMUX_DISCRETIZATION_BOX_GLOBAL_FVGEOMETRY_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/box/fvelementgeometry.hh>
#include <dumux/implicit/box/properties.hh>
#include <dumux/common/elementmap.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableGlobalFVGeometryCache>
class BoxGlobalFVGeometry
{};

// specialization in case the FVElementGeometries are stored
template<class TypeTag>
class BoxGlobalFVGeometry<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;
    //! The local class needs access to the scv, scvfs
    //! as they are globally cached
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

    using GeometryHelper = BoxGeometryHelper<GridView, dim>;

public:
    //! Constructor
    BoxGlobalFVGeometry(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView) {}

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

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.index()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! Return the gridView this global object lives on
    const GridView& gridView() const
    { return gridView_; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;

        scvs_.clear();
        scvfs_.clear();
        elementMap_.clear();

        auto numElements = gridView_.size(0);
        scvs_.resize(numElements);
        scvfs_.resize(numElements);
        elementMap_.resize(numElements);

        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        // Build the SCV and SCV faces
        for (const auto& element : elements(gridView_))
        {
            // fill the element map with seeds
            auto eIdx = problem.elementMapper().index(element);
            elementMap_[eIdx] = element.seed();

            // count
            numScv_ += element.subEntities(dim);
            numScvf_ += element.subEntities(dim-1);

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto& referenceElement = ReferenceElements::general(elementGeometry.type());

            // instantiate the geometry helper
            GeometryHelper geometryHelper(elementGeometry);

            // construct the sub control volumes
            scvs_[eIdx].reserve(elementGeometry.corners());
            for (unsigned int scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
            {
                auto dofIdxGlobal = problem.vertexMapper().subIndex(element, scvLocalIdx, dim);

                // get the corners and the volume of the scv
                auto scvCorners = geometryHelper.getScvCorners(scvLocalIdx);
                auto volume = geometryHelper.volume(scvCorners);
                scvs_[eIdx].emplace_back(std::move(scvCorners),
                                         volume,
                                         scvLocalIdx,
                                         eIdx,
                                         dofIdxGlobal);
            }

            // construct the sub control volume faces
            auto numInteriorScvfs = element.subEntities(1);
            unsigned int scvfLocalIdx = 0;
            scvfs_[eIdx].reserve(numInteriorScvfs);
            for (; scvfLocalIdx < numInteriorScvfs; ++scvfLocalIdx)
            {
                // find the global and local scv indices this scvf is belonging to
                std::vector<IndexType> localScvIndices({static_cast<IndexType>(referenceElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                        static_cast<IndexType>(referenceElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

                // get the corner points, the area, and the unit normal vector of the scv face
                auto scvfCorners = geometryHelper.getScvfCorners(scvfLocalIdx);
                auto area = geometryHelper.area(scvfCorners);
                auto normal = geometryHelper.normal(elementGeometry, scvfCorners);
                const auto v = elementGeometry.corner(localScvIndices[1]) - elementGeometry.corner(localScvIndices[0]);
                const auto s = v*normal;
                if (std::signbit(s))
                    normal *= -1;

                scvfs_[eIdx].emplace_back(std::move(scvfCorners),
                                          normal,
                                          area,
                                          scvfLocalIdx,
                                          localScvIndices,
                                          false);
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (intersection.boundary())
                {
                    auto isGeometry = intersection.geometry();
                    // count
                    numScvf_ += isGeometry.corners();
                    numBoundaryScvf_ += isGeometry.corners();

                    for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < isGeometry.corners(); ++isScvfLocalIdx)
                    {
                        // find the scvs this scvf is belonging to
                        std::vector<IndexType> localScvIndices =
                            {static_cast<IndexType>(referenceElement.subEntity(intersection.indexInInside(), 1,
                                                                               isScvfLocalIdx, dim))};

                        // get the corner points and the area of the boundary scv face
                        auto scvfCorners = geometryHelper.getBoundaryScvfCorners(isGeometry, isScvfLocalIdx);
                        auto area = geometryHelper.area(scvfCorners);

                        scvfs_[eIdx].emplace_back(std::move(scvfCorners),
                                                  intersection.centerUnitOuterNormal(),
                                                  area,
                                                  scvfLocalIdx,
                                                  localScvIndices,
                                                  true);

                        // increment local counter
                        scvfLocalIdx++;
                    }
                }
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend FVElementGeometry localView(const BoxGlobalFVGeometry& global)
    {
        return FVElementGeometry(global);
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

private:
    //! Get the local scvs for an element
    const std::vector<SubControlVolume>& scvs(IndexType eIdx) const
    { return scvs_[eIdx]; }

    //! Get the local scvfs for an element
    const std::vector<SubControlVolumeFace>& scvfs(IndexType eIdx) const
    { return scvfs_[eIdx]; }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
    const FeCache feCache_;

    GridView gridView_;
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::vector<SubControlVolume>> scvs_;
    std::vector<std::vector<SubControlVolumeFace>> scvfs_;
    // TODO do we need those?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class BoxGlobalFVGeometry<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    //! The local class needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    //! Constructor
    BoxGlobalFVGeometry(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView) {}

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

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.index()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! Return the gridView this global object lives on
    const GridView& gridView() const
    { return gridView_; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;
        elementMap_.clear();

        auto numElems = gridView_.size(0);
        elementMap_.resize(numElems);

        // save global data on the grid's scvs and scvfs
        numScv_ = 0;
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(gridView_))
        {
            // fill the element map with seeds
            auto eIdx = problem.elementMapper().index(element);
            elementMap_[eIdx] = element.seed();

            numScv_ += element.subEntities(dim);
            numScvf_ += element.subEntities(dim-1);

            // store the sub control volume face indices on the domain boundary
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (intersection.boundary())
                {
                    auto isGeometry = intersection.geometry();
                    numScvf_ += isGeometry.corners();
                    numBoundaryScvf_ += isGeometry.corners();
                }
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend FVElementGeometry localView(const BoxGlobalFVGeometry& global)
    {
        return FVElementGeometry(global);
    }

    //! The finite element cache for creating local FE bases
    const FeCache& feCache() const
    { return feCache_; }

private:

    const Problem& problem_() const
    { return *problemPtr_; }

    GridView gridView_;
    const Problem* problemPtr_;
    const FeCache feCache_;

    // Information on the global number of geometries
    // TODO do we need those information?
    std::size_t numScv_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // vectors that store the global data
    Dumux::ElementMap<GridView> elementMap_;
};

} // end namespace

#endif
