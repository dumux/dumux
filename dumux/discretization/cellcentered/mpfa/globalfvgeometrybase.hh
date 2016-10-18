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
 *        for each element.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_GLOBALFVGEOMETRY_BASE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_GLOBALFVGEOMETRY_BASE_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/common/elementmap.hh>

#include "mpfageometryhelper.hh"

namespace Dumux
{
//! Helper class for the global fv geometries
template<class TypeTag>
class CCMpfaGlobalFVGeometryHelper
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using CoordScalar = typename GridView::ctype;

    static const int dim = GridView::dimension;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    static std::vector<bool> findGhostVertices(const Problem& problem, const GridView& gridView)
    {
        std::vector<bool> ghostVertices(gridView.size(dim), false);

        // if not run in parallel, skip the rest
        if (Dune::MPIHelper::getCollectiveCommunication().size() == 1)
            return ghostVertices;

        // mpfa methods can not yet handle ghost cells
        if (gridView.ghostSize(0) > 0)
            DUNE_THROW(Dune::InvalidStateException, "Mpfa methods in parallel do not work with ghost cells. Use overlap cells instead.");

        // mpfa methods have to have overlapping cells
        if (gridView.overlapSize(0) == 0)
            DUNE_THROW(Dune::InvalidStateException, "Grid no overlaping cells. This is required by mpfa methods in parallel.");

        for (const auto& element : elements(gridView))
        {
            for (const auto& is : intersections(gridView, element))
            {
                if (!is.neighbor() && !is.boundary())
                {
                    const auto& refElement = ReferenceElements::general(element.geometry().type());
                    for (unsigned int isVertex = 0; isVertex < is.geometry().corners(); ++isVertex)
                    {
                        const auto vIdxLocal = refElement.subEntity(is.indexInInside(), 1, isVertex, dim);
                        const auto vIdxGlobal = problem.vertexMapper().subIndex(element, vIdxLocal, dim);

                        ghostVertices[vIdxGlobal] = true;
                    }
                }
            }
        }

        return ghostVertices;
    }

    template<int d = dim>
    static typename std::enable_if<d == 2, std::size_t>::type getGlobalNumScvf(const GridView& gridView)
    {
        Dune::GeometryType triangle, quadrilateral;
        triangle.makeTriangle();
        quadrilateral.makeQuadrilateral();

        return gridView.size(triangle)*6 + gridView.size(quadrilateral)*8;
    }

    template<int d = dim>
    static typename std::enable_if<d == 3, std::size_t>::type getGlobalNumScvf(const GridView& gridView)
    {
        Dune::GeometryType simplex, pyramid, prism, cube;
        simplex.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        cube.makeHexahedron();

        return gridView.size(simplex)*12 + gridView.size(pyramid)*16 + gridView.size(prism)*18 + gridView.size(cube)*24;
    }
};

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableFVElementGeometryCache>
class CCMpfaGlobalFVGeometryBase
{};

// specialization in case the FVElementGeometries are stored
template<class TypeTag>
class CCMpfaGlobalFVGeometryBase<TypeTag, true>
{
    //! The local class needs access to the scv, scvfs and the fv element geometry
    //! as they are globally cached
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using Implementation = typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GlobalInteractionVolumeSeeds = typename GET_PROP_TYPE(TypeTag, GlobalInteractionVolumeSeeds);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using BoundaryInteractionVolumeSeed = typename BoundaryInteractionVolume::Seed;
    using Element = typename GridView::template Codim<0>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;
    using IndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;
    using MpfaGeometryHelper = Dumux::MpfaGeometryHelper<GridView, dim>;

public:
    //! Constructor
    CCMpfaGlobalFVGeometryBase(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView), globalInteractionVolumeSeeds_(gridView) {}

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return scvs_.size(); }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    { return scvfs_.size(); }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! Get the inner interaction volume seed corresponding to an scvf
    const InteractionVolumeSeed& interactionVolumeSeed(const SubControlVolumeFace& scvf) const
    { return globalInteractionVolumeSeeds_.seed(scvf); }

    //! Get the boundary interaction volume seed corresponding to an scvf
    const BoundaryInteractionVolumeSeed& boundaryInteractionVolumeSeed(const SubControlVolumeFace& scvf) const
    { return globalInteractionVolumeSeeds_.boundarySeed(scvf); }

    //! Returns whether or not a scvf touches the boundary (has to be called before getting an interaction volume)
    bool scvfTouchesBoundary(const SubControlVolumeFace& scvf) const
    { return boundaryVertices_[scvf.vertexIndex()]; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;

        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvfIndicesOfScv_.clear();
        boundaryVertices_.clear();
        elementMap_.clear();

        // reserve memory
        std::size_t numScvs = gridView_.size(0);
        std::size_t numScvfs = CCMpfaGlobalFVGeometryHelper<TypeTag>::getGlobalNumScvf(gridView_);
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvfs);
        scvfIndicesOfScv_.resize(numScvs);
        boundaryVertices_.resize(gridView_.size(dim), false);
        elementMap_.resize(numScvs);

        // find vertices on processor boundaries
        auto ghostVertices = CCMpfaGlobalFVGeometryHelper<TypeTag>::findGhostVertices(problem, gridView_);

        // the quadrature point to be used on the scvf
        const Scalar q = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, Q);

        // Build the SCVs and SCV faces
        IndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(gridView_))
        {
            auto eIdx = problem.elementMapper().index(element);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element geometry and reference element
            auto elemGeometry = element.geometry();
            const auto& referenceElement = ReferenceElements::general(elemGeometry.type());

            // The geometry helper class
            MpfaGeometryHelper geomHelper(elemGeometry);

            // create sub control volume for this element
            scvs_[eIdx] = SubControlVolume(std::move(elemGeometry), eIdx);

            // The local scvf index set
            std::vector<IndexType> scvfIndexSet;
            scvfIndexSet.reserve(geomHelper.getNumLocalScvfs());

            // construct the sub control volume faces
            for (const auto& intersection : intersections(gridView_, element))
            {
                // get the intersection geometry and some of its bools
                auto isGeometry = intersection.geometry();
                bool boundary = intersection.boundary();

                // determine the outside volvar idx
                IndexType nIdx;
                if (intersection.neighbor())
                    nIdx = problem.elementMapper().index(intersection.outside());
                else if (boundary)
                    nIdx = numScvs + numBoundaryScvf_++;

                // make the scv faces of the intersection
                for (unsigned int faceScvfIdx = 0; faceScvfIdx < isGeometry.corners(); ++faceScvfIdx)
                {
                    // get the global vertex index the scv face is connected to (mpfa-o method does not work for hanging nodes!)
                    const auto vIdxLocal = referenceElement.subEntity(intersection.indexInInside(), 1, faceScvfIdx, dim);
                    const auto vIdxGlobal = problem.vertexMapper().subIndex(element, vIdxLocal, dim);

                    // do not build scvfs connected to a processor boundary
                    if (ghostVertices[vIdxGlobal])
                        continue;

                    // store info on which vertices are on the domain boundary
                    if (boundary)
                        boundaryVertices_[vIdxGlobal] = true;

                    // make the scv face
                    scvfIndexSet.push_back(scvfIdx);
                    scvfs_.emplace_back(geomHelper,
                                        geomHelper.getScvfCorners(isGeometry, faceScvfIdx),
                                        intersection.centerUnitOuterNormal(),
                                        vIdxGlobal,
                                        scvfIdx,
                                        std::array<IndexType, 2>({{eIdx, nIdx}}),
                                        q,
                                        boundary
                                        );

                    // increment scvf counter
                    scvfIdx++;
                }
            }

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfIndexSet;
        }

        // the number of actual boundary scvf is two times the number of boundary intersections
        numBoundaryScvf_ *= 2;

        // in parallel problems we might have reserved more scvfs than we actually use
        scvfs_.shrink_to_fit();

        // Initialize the interaction volume seeds
        globalInteractionVolumeSeeds_.update(problem, boundaryVertices_);
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline FVElementGeometry localView(const Implementation& global)
    { return FVElementGeometry(global); }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(IndexType scvIdx) const
    { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

private:

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
    GridView gridView_;

    // vectors that store the geometries
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;

    // vectors that store the global data
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<bool> boundaryVertices_;
    IndexType numBoundaryScvf_;

    // the global interaction volume seeds
    GlobalInteractionVolumeSeeds globalInteractionVolumeSeeds_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class CCMpfaGlobalFVGeometryBase<TypeTag, false>
{
    //! The local fvGeometry needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using Implementation = typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GlobalInteractionVolumeSeeds = typename GET_PROP_TYPE(TypeTag, GlobalInteractionVolumeSeeds);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using BoundaryInteractionVolumeSeed = typename BoundaryInteractionVolume::Seed;
    using Element = typename GridView::template Codim<0>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;
    using IndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;
    using MpfaGeometryHelper = Dumux::MpfaGeometryHelper<GridView, dim>;

public:
    //! Constructor
    CCMpfaGlobalFVGeometryBase(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView), globalInteractionVolumeSeeds_(gridView_) {}

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return numScvs_; }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! Return the gridView this global object lives on
    const GridView& gridView() const
    { return gridView_; }

    //! Get the inner interaction volume seed corresponding to an scvf
    const InteractionVolumeSeed& interactionVolumeSeed(const SubControlVolumeFace& scvf) const
    { return globalInteractionVolumeSeeds_.seed(scvf); }

    //! Get the boundary interaction volume seed corresponding to an scvf
    const BoundaryInteractionVolumeSeed& boundaryInteractionVolumeSeed(const SubControlVolumeFace& scvf) const
    { return globalInteractionVolumeSeeds_.boundarySeed(scvf); }

    //! Returns whether or not a scvf touches the boundary (has to be called before getting an interaction volume)
    bool scvfTouchesBoundary(const SubControlVolumeFace& scvf) const
    { return boundaryVertices_[scvf.vertexIndex()]; }

    //! Returns whether or not a vertex is on a processor boundary
    bool isGhostVertex(const IndexType vIdxGlobal) const
    { return ghostVertices_[vIdxGlobal]; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;

        // clear containers (necessary after grid refinement)
        scvfIndicesOfScv_.clear();
        neighborVolVarIndices_.clear();
        boundaryVertices_.clear();
        elementMap_.clear();

        // reserve memory or resize the containers
        numScvs_ = gridView_.size(0);
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        elementMap_.resize(numScvs_);
        scvfIndicesOfScv_.resize(numScvs_);
        neighborVolVarIndices_.resize(numScvs_);
        boundaryVertices_.resize(gridView_.size(dim), false);

        // find vertices on processor boundaries
        ghostVertices_ = CCMpfaGlobalFVGeometryHelper<TypeTag>::findGhostVertices(problem, gridView_);

        // Store necessary info on SCVs and SCV faces
        for (const auto& element : elements(gridView_))
        {
            auto eIdx = problem.elementMapper().index(element);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element geometry and reference element
            auto elemGeometry = element.geometry();
            const auto& referenceElement = ReferenceElements::general(elemGeometry.type());

            // The geometry helper class
            MpfaGeometryHelper geomHelper(elemGeometry);

            // the element-wise index sets for finite volume geometry
            auto numLocalFaces = geomHelper.getNumLocalScvfs();
            std::vector<IndexType> scvfsIndexSet(numLocalFaces);
            std::vector<IndexType> neighborVolVarIndexSet(numLocalFaces);
            IndexType localFaceIdx = 0;
            // construct the sub control volume faces
            for (const auto& intersection : intersections(gridView_, element))
            {
                // get some of the intersection bools
                bool boundary = intersection.boundary();
                bool neighbor = intersection.neighbor();

                // determine the outside volvar idx
                IndexType nIdx;
                if (neighbor)
                    nIdx = problem.elementMapper().index(intersection.outside());
                else if (boundary)
                    nIdx = numScvs_ + numBoundaryScvf_++;

                // make the scv faces of the intersection
                for (unsigned int faceScvfIdx = 0; faceScvfIdx < intersection.geometry().corners(); ++faceScvfIdx)
                {
                    // get the global vertex index the scv face is connected to (mpfa-o method does not work for hanging nodes!)
                    const auto vIdxLocal = referenceElement.subEntity(intersection.indexInInside(), 1, faceScvfIdx, dim);
                    const auto vIdxGlobal = problem.vertexMapper().subIndex(element, vIdxLocal, dim);

                    // do not build scvfs connected to a processor boundary
                    if (ghostVertices_[vIdxGlobal])
                        continue;

                    // store info on which vertices are on the domain boundary
                    if (boundary)
                        boundaryVertices_[vIdxGlobal] = true;

                    // store information on the scv face
                    scvfsIndexSet[localFaceIdx] = numScvf_++;
                    neighborVolVarIndexSet[localFaceIdx++] = nIdx;
                }
            }

            // store the sets of indices in the data container
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
            neighborVolVarIndices_[eIdx] = neighborVolVarIndexSet;
        }

        // the number of actual boundary scvf is two times the number of boundary intersections
        numBoundaryScvf_ *= 2;

        // Initialize the interaction volume seeds
        globalInteractionVolumeSeeds_.update(problem, boundaryVertices_);
    }

    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    const std::vector<IndexType>& neighborVolVarIndices(IndexType scvIdx) const
    { return neighborVolVarIndices_[scvIdx]; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline FVElementGeometry localView(const Implementation& global)
    { return FVElementGeometry(global); }

private:

    //! Get the inner interaction volume seed corresponding to an scvf
    const InteractionVolumeSeed& interactionVolumeSeed(const IndexType vIdxGlobal) const
    { return globalInteractionVolumeSeeds_.seed(vIdxGlobal); }

    //! Get the boundary interaction volume seed corresponding to an scvf
    const BoundaryInteractionVolumeSeed& boundaryInteractionVolumeSeed(const IndexType vIdxGlobal) const
    { return globalInteractionVolumeSeeds_.boundarySeed(vIdxGlobal); }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    GridView gridView_;

    // Information on the global number of geometries
    IndexType numScvs_;
    IndexType numScvf_;
    IndexType numBoundaryScvf_;

    // vectors that store the global data
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<IndexType>> neighborVolVarIndices_;
    std::vector<bool> boundaryVertices_;
    std::vector<bool> ghostVertices_;

    // the global interaction volume seeds
    GlobalInteractionVolumeSeeds globalInteractionVolumeSeeds_;
};

} // end namespace

#endif
