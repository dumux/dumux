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
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_MIMETIC_GLOBAL_FVGEOMETRY_HH
#define DUMUX_DISCRETIZATION_MIMETIC_GLOBAL_FVGEOMETRY_HH

#include <dumux/common/elementmap.hh>
#include <dumux/implicit/staggered/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableGlobalFVGeometryCache>
class MimeticGlobalFVGeometry
{};

// specialization in case the FVElementGeometries are stored globally
template<class TypeTag>
class MimeticGlobalFVGeometry<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;
    using IntersectionMapper = typename GET_PROP_TYPE(TypeTag, IntersectionMapper);
    //! The local class needs access to the scv, scvfs and the fv element geometry
    //! as they are globally cached
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using GeometryHelper = typename GET_PROP_TYPE(TypeTag, StaggeredGeometryHelper);

public:
    //! Constructor
    MimeticGlobalFVGeometry(const GridView& gridView)
    : gridView_(gridView), elementMap_(gridView), intersectionMapper_(gridView) {}

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
    }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    {
        return numBoundaryScvf_;
    }

    //! The total number of intersections
    std::size_t numIntersections() const
    {
        return intersectionMapper_.numIntersections();
    }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

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

        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvfIndicesOfScv_.clear();
        elementMap_.clear();
        intersectionMapper_.update();

        // determine size of containers
        IndexType numScvs = gridView_.size(0);
        IndexType numScvf = 0;
        for (const auto& element : elements(gridView_))
            numScvf += element.subEntities(1);

        // reserve memory
        elementMap_.resize(numScvs);
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvf);
        scvfIndicesOfScv_.resize(numScvs);
        localToGlobalScvfIndices_.resize(numScvs);

        // Build the scvs and scv faces
        IndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(gridView_))
        {
            bool calcNewCellCenter = false;
            auto eIdx = problem.elementMapper().index(element);

            // reserve memory for the localToGlobalScvfIdx map
            auto numLocalFaces = intersectionMapper_.numFaces(element);
            localToGlobalScvfIndices_[eIdx].resize(numLocalFaces);

            scvs_[eIdx] = SubControlVolume(element.geometry(), eIdx);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);

            GeometryHelper geometryHelper(element, gridView_);

            for (const auto& intersection : intersections(gridView_, element))
            {
                geometryHelper.updateLocalFace(intersectionMapper_, intersection);
                const int localFaceIndex = geometryHelper.localFaceIndex();

                // inner sub control volume faces
                if (intersection.neighbor())
                {
                    auto nIdx = problem.elementMapper().index(intersection.outside());
                    scvfs_.emplace_back(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<IndexType>({eIdx, nIdx}),
                                        geometryHelper
                                        );
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;

                    auto di = scvfs_[scvfIdx].ipGlobal();
                    di -= element.geometry().center();
                    if(scvfs_[scvfIdx].unitOuterNormal()*di < 0)
                        calcNewCellCenter = true;

                    scvfsIndexSet.push_back(scvfIdx++);
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    scvfs_.emplace_back(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<IndexType>({eIdx, gridView_.size(0) + numBoundaryScvf_++}),
                                        geometryHelper
                                        );
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;

                    auto di = scvfs_[scvfIdx].ipGlobal();
                    di -= element.geometry().center();
                    if(scvfs_[scvfIdx].unitOuterNormal()*di < 0)
                        calcNewCellCenter = true;

                    scvfsIndexSet.push_back(scvfIdx++);
                }
            }

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;

            //if(calcNewCellCenter)
            //    findNewCellCenter(eIdx);
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline FVElementGeometry localView(const MimeticGlobalFVGeometry& global)
    { return FVElementGeometry(global); }

//private:

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        return scvfs_[scvfIdx];
    }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    {
        return scvfIndicesOfScv_[scvIdx];
    }

    const auto localToGlobalScvfIndex(IndexType eIdx, IndexType localScvfIdx) const
    {
        return localToGlobalScvfIndices_[eIdx][localScvfIdx];
    }

    const SubControlVolumeFace& scvf(IndexType eIdx ,IndexType localScvfIdx) const
    {
        return scvf(localToGlobalScvfIndex(eIdx, localScvfIdx));
    }

    void findNewCellCenter(int eIdx)
    {
        GlobalPosition center(0);
        int fIdxI = -1;
        int fIdxJ = -1;
        double distFuncVal = 1.0e100;
        int numFaces = localToGlobalScvfIndices_[eIdx].size();

        for(int i=0; i<numFaces; i++)
        {
            auto scvfIdx = localToGlobalScvfIndices_[eIdx][i];
            const auto scfv = scvfs_[scvfIdx];
            auto faceCenterI = scfv.ipGlobal();
            for(int j=0; j<numFaces; j++)
            {
                auto scvfIdxJ = localToGlobalScvfIndices_[eIdx][j];
                const auto scfvJ = scvfs_[scvfIdxJ];
                GlobalPosition faceCenterJ = scfvJ.ipGlobal();
                center = faceCenterI + faceCenterJ;
                center /= 2.0;
                if(checkValidility(center, eIdx))
                {
                    Scalar distVal = calculateInvDistFunction(center, eIdx);
                    if(distVal < distFuncVal)
                    {
                        fIdxI = i;
                        fIdxJ = j;
                        distFuncVal = distVal;
                    }
                }
            }

        }

        if(fIdxI == -1 && fIdxJ == -1 )
            DUNE_THROW(Dune::InvalidStateException, "Cannot find new cell center!");
        else
        {
            auto scvfIdx = localToGlobalScvfIndices_[eIdx][fIdxI];
            const auto scfvI = scvfs_[scvfIdx];
            auto scvfIdxJ = localToGlobalScvfIndices_[eIdx][fIdxJ];
            const auto scfvJ = scvfs_[scvfIdxJ];
            GlobalPosition faceCenterI = scfvI.ipGlobal();
            GlobalPosition faceCenterJ = scfvJ.ipGlobal();
            center = faceCenterI + faceCenterJ;
            center /= 2.0;
            std::cout << "Found new cell center: " << center << " Old center was: " << scvs_[eIdx].center() << std::endl;

            scvs_[eIdx].setCellCenter(center);
        }

    }

    bool checkValidility(GlobalPosition& Point, int eIdx)
    {
        for(int i=0; i<localToGlobalScvfIndices_[eIdx].size(); i++)
        {
            auto scvfIdx = localToGlobalScvfIndices_[eIdx][i];
            const auto scfv = scvfs_[scvfIdx];
            if(scfv.unitOuterNormal() * (scfv.ipGlobal()- Point) < 1.0e-8)
                return false;
        }

        return true;
    }

    double calculateInvDistFunction(GlobalPosition& Point, int eIdx)
    {
        double distFuncVal = 0.0;
        for(int i=0; i<localToGlobalScvfIndices_[eIdx].size(); i++)
        {
            auto scvfIdx = localToGlobalScvfIndices_[eIdx][i];
            const auto scfv = scvfs_[scvfIdx];
            Scalar val = scfv.unitOuterNormal() * (scfv.ipGlobal()- Point);
            distFuncVal += 1.0/val*val;
        }
        return distFuncVal;
    }


public:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    const GridView gridView_;
    Dumux::ElementMap<GridView> elementMap_;
    IntersectionMapper intersectionMapper_;
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<IndexType>> localToGlobalScvfIndices_;
    IndexType numBoundaryScvf_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class MimeticGlobalFVGeometry<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;
    //! The local fvGeometry needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

public:
    //! Constructor
    MimeticGlobalFVGeometry(const GridView& gridView)
    : gridView_(gridView), elementMap_(gridView) {}

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return numScvs_;
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return numScvf_;
    }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    {
        return numBoundaryScvf_;
    }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

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

        // reserve memory or resize the containers
        numScvs_ = gridView_.size(0);
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        elementMap_.resize(numScvs_);
        scvfIndicesOfScv_.resize(numScvs_);
        neighborVolVarIndices_.resize(numScvs_);

        // Build the SCV and SCV face
        for (const auto& element : elements(gridView_))
        {
            auto eIdx = problem.elementMapper().index(element);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element-wise index sets for finite volume geometry
            auto numLocalFaces = element.subEntities(1);
            std::vector<IndexType> scvfsIndexSet;
            std::vector<IndexType> neighborVolVarIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);
            neighborVolVarIndexSet.reserve(numLocalFaces);
            for (const auto& intersection : intersections(gridView_, element))
            {
                // inner sub control volume faces
                if (intersection.neighbor())
                {
                    scvfsIndexSet.push_back(numScvf_++);
                    neighborVolVarIndexSet.push_back(problem.elementMapper().index(intersection.outside()));
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    scvfsIndexSet.push_back(numScvf_++);
                    neighborVolVarIndexSet.push_back(numScvs_ + numBoundaryScvf_++);
                }
            }

            // store the sets of indices in the data container
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
            neighborVolVarIndices_[eIdx] = neighborVolVarIndexSet;
        }
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
    friend inline FVElementGeometry localView(const MimeticGlobalFVGeometry& global)
    { return FVElementGeometry(global); }

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    const GridView gridView_;

    // Information on the global number of geometries
    IndexType numScvs_;
    IndexType numScvf_;
    IndexType numBoundaryScvf_;

    // vectors that store the global data
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<IndexType>> neighborVolVarIndices_;
};

} // end namespace

#endif
