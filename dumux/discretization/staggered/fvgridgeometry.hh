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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FV_GRID_GEOMETRY
#define DUMUX_DISCRETIZATION_STAGGERED_FV_GRID_GEOMETRY

#include <dumux/discretization/basefvgridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for cell center of face specific auxiliary FvGridGeometry classes.
 *        Provides a common interface and a pointer to the actual fvGridGeometry.
 */
template<class ActualFVGridGeometry>
class GridGeometryView
{
public:

    explicit GridGeometryView(const ActualFVGridGeometry* actualFVGridGeometry)
    : fvGridGeometry_(actualFVGridGeometry) {}

    //! export  the GridView type and the discretization method
    using GridView = typename ActualFVGridGeometry::GridView;
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::staggered;
    using LocalView = typename ActualFVGridGeometry::LocalView;

    /*!
     * \brief Returns true if this view if related to cell centered dofs
     */
    static constexpr bool isCellCenter() { return false; }

    /*!
     * \brief Returns true if this view if related to face dofs
     */
    static constexpr bool isFace() {return false; }

    /*!
     * \brief Return an integral constant index for cell centered dofs
     */
    static constexpr auto cellCenterIdx()
    { return typename ActualFVGridGeometry::DofTypeIndices::CellCenterIdx{}; }

    /*!
     * \brief Return an integral constant index for face dofs
     */
    static constexpr auto faceIdx()
    { return typename ActualFVGridGeometry::DofTypeIndices::FaceIdx{}; }

    /*!
     * \brief Return the gridView this grid geometry object lives on
     */
    const auto& gridView() const
    { return fvGridGeometry_->gridView(); }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const auto& connectivityMap() const // TODO return correct map
    { return fvGridGeometry_->connectivityMap(); }

    /*!
     * \brief Returns the mapper for vertices to indices for possibly adaptive grids.
     */
    const auto& vertexMapper() const
    { return fvGridGeometry_->vertexMapper(); }

    /*!
     * \brief Returns the mapper for elements to indices for constant grids.
     */
    const auto& elementMapper() const
    { return fvGridGeometry_->elementMapper(); }

    /*!
     * \brief Returns the actual fvGridGeometry we are a restriction of
     */
    const ActualFVGridGeometry& actualfvGridGeometry() const
    { return *fvGridGeometry_; }

protected:
    const ActualFVGridGeometry* fvGridGeometry_;

};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Cell center specific auxiliary FvGridGeometry classes.
 *        Required for the Dumux multi-domain framework.
 */
template <class ActualFVGridGeometry>
class CellCenterFVGridGeometry : public GridGeometryView<ActualFVGridGeometry>
{
    using ParentType = GridGeometryView<ActualFVGridGeometry>;
public:

    using ParentType::ParentType;

    /*!
     * \brief Returns true because this view is related to cell centered dofs
     */
    static constexpr bool isCellCenter() { return true; }

    /*!
     * \brief The total number of cell centered dofs
     */
    std::size_t numDofs() const
    { return this->fvGridGeometry_->numCellCenterDofs(); }
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face specific auxiliary FvGridGeometry classes.
 *        Required for the Dumux multi-domain framework.
 */
template <class ActualFVGridGeometry>
class FaceFVGridGeometry : public GridGeometryView<ActualFVGridGeometry>
{
    using ParentType = GridGeometryView<ActualFVGridGeometry>;
public:

    using ParentType::ParentType;

    /*!
     * \brief Returns true because this view is related to face dofs
     */
    static constexpr bool isFace() {return true; }

    /*!
     * \brief The total number of cell centered dofs
     */
    std::size_t numDofs() const
    { return this->fvGridGeometry_->numFaceDofs(); }
};



/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class GridView,
         bool cachingEnabled,
         class Traits>
class StaggeredFVGridGeometry;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored.
 */
template<class GV, class Traits>
class StaggeredFVGridGeometry<GV, true, Traits>
: public BaseFVGridGeometry<StaggeredFVGridGeometry<GV, true, Traits>, GV, Traits>
{
    using ThisType = StaggeredFVGridGeometry<GV, true, Traits>;
    using ParentType = BaseFVGridGeometry<ThisType, GV, Traits>;
    using IndexType = typename GV::IndexSet::IndexType;

    using IntersectionMapper = typename Traits::IntersectionMapper;
    using GeometryHelper = typename Traits::GeometryHelper;
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::staggered;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the grid view type
    using GridView = GV;
    using Element = typename GV::template Codim<0>::Entity;
    //! export the dof type indices
    using DofTypeIndices = typename Traits::DofTypeIndices;

    //! return a integral constant for cell center dofs
    static constexpr auto cellCenterIdx()
    { return typename DofTypeIndices::CellCenterIdx{}; }

    //! return a integral constant for face dofs
    static constexpr auto faceIdx()
    { return typename DofTypeIndices::FaceIdx{}; }

    using CellCenterFVGridGeometryType = CellCenterFVGridGeometry<ThisType>;
    using FaceFVGridGeometryType = FaceFVGridGeometry<ThisType>;

    using FVGridGeometryTuple = std::tuple< CellCenterFVGridGeometry<ThisType>, FaceFVGridGeometry<ThisType> >;

    //! Constructor
    StaggeredFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    , intersectionMapper_(gridView)
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod::staggered>::isValid(gridView))
            DUNE_THROW(Dune::InvalidStateException, "The staggered discretization method needs at least an overlap of 1 for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");
    }

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

    //! the total number of dofs
    std::size_t numDofs() const
    { return numCellCenterDofs() + numFaceDofs(); }

    std::size_t numCellCenterDofs() const
    { return this->gridView().size(0); }

    std::size_t numFaceDofs() const
    { return this->gridView().size(1); }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update()
    {
        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        boundaryScvfs_.clear();
        scvfIndicesOfScv_.clear();
        intersectionMapper_.update();

        // determine size of containers
        IndexType numScvs = this->gridView().size(0);
        IndexType numScvf = 0;
        for (const auto& element : elements(this->gridView()))
            numScvf += element.subEntities(1);

        // reserve memory
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvf);
        boundaryScvfs_.resize(numScvf);
        scvfIndicesOfScv_.resize(numScvs);
        localToGlobalScvfIndices_.resize(numScvs);

        // Build the scvs and scv faces
        IndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            auto eIdx = this->elementMapper().index(element);

            // reserve memory for the localToGlobalScvfIdx map
            auto numLocalFaces = intersectionMapper_.numFaces(element);
            localToGlobalScvfIndices_[eIdx].resize(numLocalFaces);

            scvs_[eIdx] = SubControlVolume(element.geometry(), eIdx);

            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);

            GeometryHelper geometryHelper(element, this->gridView());

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                geometryHelper.updateLocalFace(intersectionMapper_, intersection);
                const int localFaceIndex = geometryHelper.localFaceIndex();

                // inner sub control volume faces
                if (intersection.neighbor())
                {
                    auto nIdx = this->elementMapper().index(intersection.outside());
                    scvfs_.emplace_back(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<IndexType>({eIdx, nIdx}),
                                        geometryHelper);
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;
                    scvfsIndexSet.push_back(scvfIdx++);
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    SubControlVolumeFace tmpScvf(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<IndexType>({eIdx, this->gridView().size(0) + numBoundaryScvf_}),
                                        geometryHelper
                                        );
                    boundaryScvfs_[this->gridView().indexSet().subIndex(element, localFaceIndex, 1)] = tmpScvf;
                    scvfs_.emplace_back(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<IndexType>({eIdx, this->gridView().size(0) + numBoundaryScvf_++}),
                                        geometryHelper);
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;
                    scvfsIndexSet.push_back(scvfIdx++);
                    boundaryScvfsIndexSet_.push_back(this->gridView().indexSet().subIndex(element, localFaceIndex, 1));
                    bool eIdxAlreadyInVector = (find(boundaryScvsIndexSet_.begin(), boundaryScvsIndexSet_.end(), eIdx) != boundaryScvsIndexSet_.end());
                    if(eIdxAlreadyInVector == false)
                        boundaryScvsIndexSet_.push_back(eIdx);
                }
            }

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
        }

        // build the connectivity map for an effecient assembly
        connectivityMap_.update(*this);
    }

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

    //! Get a boundary sub control volume face with dof index
    // At the boundary, the dofIdx is clearly connected with a scvf. In the inner of the domain, the dofIdx belongs to two different scvfs. ScvfIdx is not dofIdx, scvfIdx is different depending on what element the scvf is associated with.
    const SubControlVolumeFace& boundaryScvf(IndexType scvfDofIdx) const
    {
        return boundaryScvfs_[scvfDofIdx];
    }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    {
        return scvfIndicesOfScv_[scvIdx];
    }

    //! Get a vector of all dofIndices of the boundary scvfs
    const std::vector<IndexType>& boundaryScvfsIndexSet() const
    {
        return boundaryScvfsIndexSet_;
    }

    //! Get a vector of all dofIndices of the boundary scvs
    const std::vector<IndexType>& boundaryScvsIndexSet() const
    {
        return boundaryScvsIndexSet_;
    }

    IndexType localToGlobalScvfIndex(IndexType eIdx, IndexType localScvfIdx) const
    {
        return localToGlobalScvfIndices_[eIdx][localScvfIdx];
    }

    const SubControlVolumeFace& scvf(IndexType eIdx ,IndexType localScvfIdx) const
    {
        return scvf(localToGlobalScvfIndex(eIdx, localScvfIdx));
    }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const ConnectivityMap &connectivityMap() const
    { return connectivityMap_; }

    //! Returns a pointer the cell center specific auxiliary class. Required for the multi-domain FVAssembler's ctor.
    std::unique_ptr<CellCenterFVGridGeometry<ThisType>> cellCenterFVGridGeometryPtr() const
    {
        return std::make_unique<CellCenterFVGridGeometry<ThisType>>(this);
    }

    //! Returns a pointer the face specific auxiliary class. Required for the multi-domain FVAssembler's ctor.
    std::unique_ptr<FaceFVGridGeometry<ThisType>> faceFVGridGeometryPtr() const
    {
        return std::make_unique<FaceFVGridGeometry<ThisType>>(this);
    }

    //! Return a copy of the cell center specific auxiliary class.
    CellCenterFVGridGeometry<ThisType> cellCenterFVGridGeometry() const
    {
        return CellCenterFVGridGeometry<ThisType>(this);
    }

    //! Return a copy of the face specific auxiliary class.
    FaceFVGridGeometry<ThisType> faceFVGridGeometry() const
    {
        return FaceFVGridGeometry<ThisType>(this);
    }

private:

    // mappers
    ConnectivityMap connectivityMap_;
    IntersectionMapper intersectionMapper_;

    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<SubControlVolumeFace> boundaryScvfs_;
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<IndexType>> localToGlobalScvfIndices_;
    std::vector<IndexType> boundaryScvfsIndexSet_;
    std::vector<IndexType> boundaryScvsIndexSet_;
    IndexType numBoundaryScvf_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored.
 */
template<class GV, class Traits>
class StaggeredFVGridGeometry<GV, false, Traits>
: public BaseFVGridGeometry<StaggeredFVGridGeometry<GV, false, Traits>, GV, Traits>
{
    using ThisType = StaggeredFVGridGeometry<GV, false, Traits>;
    using ParentType = BaseFVGridGeometry<ThisType, GV, Traits>;
    using IndexType = typename GV::IndexSet::IndexType;

    using IntersectionMapper = typename Traits::IntersectionMapper;
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::staggered;

    using GeometryHelper = typename Traits::GeometryHelper;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, false>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the grid view type
    using GridView = GV;
    using Element = typename GV::template Codim<0>::Entity;
    //! export the dof type indices
    using DofTypeIndices = typename Traits::DofTypeIndices;

    //! return a integral constant for cell center dofs
    static constexpr auto cellCenterIdx()
    { return typename DofTypeIndices::CellCenterIdx{}; }

    //! return a integral constant for face dofs
    static constexpr auto faceIdx()
    { return typename DofTypeIndices::FaceIdx{}; }

    using CellCenterFVGridGeometryType = CellCenterFVGridGeometry<ThisType>;
    using FaceFVGridGeometryType = FaceFVGridGeometry<ThisType>;

    using FVGridGeometryTuple = std::tuple< CellCenterFVGridGeometry<ThisType>, FaceFVGridGeometry<ThisType> >;

    //! Constructor
    StaggeredFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    , intersectionMapper_(gridView)
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod::staggered>::isValid(gridView))
            DUNE_THROW(Dune::InvalidStateException, "The staggered discretization method needs at least an overlap of 1 for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");
    }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update()
    {
        // clear containers (necessary after grid refinement)
        scvfIndicesOfScv_.clear();
        intersectionMapper_.update();
        neighborVolVarIndices_.clear();

        numScvs_ = numCellCenterDofs();
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        scvfIndicesOfScv_.resize(numScvs_);
        localToGlobalScvfIndices_.resize(numScvs_);
        neighborVolVarIndices_.resize(numScvs_);

        // Build the scvs and scv faces
        for (const auto& element : elements(this->gridView()))
        {
            auto eIdx = this->elementMapper().index(element);

            // the element-wise index sets for finite volume geometry
            auto numLocalFaces = intersectionMapper_.numFaces(element);
            std::vector<IndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);
            localToGlobalScvfIndices_[eIdx].resize(numLocalFaces);

            std::vector<IndexType> neighborVolVarIndexSet;
            neighborVolVarIndexSet.reserve(numLocalFaces);

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                const auto localFaceIndex = intersection.indexInInside();
                localToGlobalScvfIndices_[eIdx][localFaceIndex] = numScvf_;
                scvfsIndexSet.push_back(numScvf_++);

                if (intersection.neighbor())
                {
                    const auto nIdx = this->elementMapper().index(intersection.outside());
                    neighborVolVarIndexSet.emplace_back(nIdx);
                }
                else
                    neighborVolVarIndexSet.emplace_back(numScvs_ + numBoundaryScvf_++);
            }

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
            neighborVolVarIndices_[eIdx] = neighborVolVarIndexSet;
        }

        // build the connectivity map for an effecient assembly
        connectivityMap_.update(*this);
    }

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

    //! The total number of intersections
    std::size_t numIntersections() const
    {
        return intersectionMapper_.numIntersections();
    }

    //! the total number of dofs
    std::size_t numDofs() const
    { return numCellCenterDofs() + numFaceDofs(); }

    std::size_t numCellCenterDofs() const
    { return this->gridView().size(0); }

    std::size_t numFaceDofs() const
    { return this->gridView().size(1); }

    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    IndexType localToGlobalScvfIndex(IndexType eIdx, IndexType localScvfIdx) const
    {
        return localToGlobalScvfIndices_[eIdx][localScvfIdx];
    }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const ConnectivityMap &connectivityMap() const
    { return connectivityMap_; }

    //! Returns a pointer the cell center specific auxiliary class. Required for the multi-domain FVAssembler's ctor.
    std::unique_ptr<CellCenterFVGridGeometry<ThisType>> cellCenterFVGridGeometryPtr() const
    {
        return std::make_unique<CellCenterFVGridGeometry<ThisType>>(this);
    }

    //! Returns a pointer the face specific auxiliary class. Required for the multi-domain FVAssembler's ctor.
    std::unique_ptr<FaceFVGridGeometry<ThisType>> faceFVGridGeometryPtr() const
    {
        return std::make_unique<FaceFVGridGeometry<ThisType>>(this);
    }

    //! Return a copy of the cell center specific auxiliary class.
    CellCenterFVGridGeometry<ThisType> cellCenterFVGridGeometry() const
    {
        return CellCenterFVGridGeometry<ThisType>(this);
    }

    //! Return a copy of the face specific auxiliary class.
    FaceFVGridGeometry<ThisType> faceFVGridGeometry() const
    {
        return FaceFVGridGeometry<ThisType>(this);
    }

    //! Return a reference to the intersection mapper
    const IntersectionMapper& intersectionMapper() const
    {
        return intersectionMapper_;
    }

    //! Return the neighbor volVar indices for all scvfs in the scv with index scvIdx
    const std::vector<IndexType>& neighborVolVarIndices(IndexType scvIdx) const
    { return neighborVolVarIndices_[scvIdx]; }

private:

    //! Information on the global number of geometries
    IndexType numScvs_;
    IndexType numScvf_;
    IndexType numBoundaryScvf_;
    std::vector<std::vector<IndexType>> localToGlobalScvfIndices_;
    std::vector<std::vector<IndexType>> neighborVolVarIndices_;

    // mappers
    ConnectivityMap connectivityMap_;
    IntersectionMapper intersectionMapper_;

    //! vectors that store the global data
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
};

} // end namespace

#endif
