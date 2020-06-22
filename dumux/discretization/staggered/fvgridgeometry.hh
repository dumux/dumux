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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FV_GRID_GEOMETRY
#define DUMUX_DISCRETIZATION_STAGGERED_FV_GRID_GEOMETRY

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for cell center of face specific auxiliary FvGridGeometry classes.
 *        Provides a common interface and a pointer to the actual gridGeometry.
 */
template<class ActualGridGeometry>
class GridGeometryView
{
public:

    explicit GridGeometryView(const ActualGridGeometry* actualGridGeometry)
    : gridGeometry_(actualGridGeometry) {}

    //! export  the GridView type and the discretization method
    using GridView = typename ActualGridGeometry::GridView;
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::staggered;
    using LocalView = typename ActualGridGeometry::LocalView;

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
    { return typename ActualGridGeometry::DofTypeIndices::CellCenterIdx{}; }

    /*!
     * \brief Return an integral constant index for face dofs
     */
    static constexpr auto faceIdx()
    { return typename ActualGridGeometry::DofTypeIndices::FaceIdx{}; }

    /*!
     * \brief Return the gridView this grid geometry object lives on
     */
    const auto& gridView() const
    { return gridGeometry_->gridView(); }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const auto& connectivityMap() const // TODO return correct map
    { return gridGeometry_->connectivityMap(); }

    /*!
     * \brief Returns the mapper for vertices to indices for possibly adaptive grids.
     */
    const auto& vertexMapper() const
    { return gridGeometry_->vertexMapper(); }

    /*!
     * \brief Returns the mapper for elements to indices for constant grids.
     */
    const auto& elementMapper() const
    { return gridGeometry_->elementMapper(); }

    /*!
     * \brief Returns the actual gridGeometry we are a restriction of
     */
    const ActualGridGeometry& actualGridGeometry() const
    { return *gridGeometry_; }

protected:
    const ActualGridGeometry* gridGeometry_;

};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Cell center specific auxiliary FvGridGeometry classes.
 *        Required for the Dumux multi-domain framework.
 */
template <class ActualGridGeometry>
class CellCenterFVGridGeometry : public GridGeometryView<ActualGridGeometry>
{
    using ParentType = GridGeometryView<ActualGridGeometry>;
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
    { return this->gridGeometry_->numCellCenterDofs(); }
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Face specific auxiliary FvGridGeometry classes.
 *        Required for the Dumux multi-domain framework.
 */
template <class ActualGridGeometry>
class FaceFVGridGeometry : public GridGeometryView<ActualGridGeometry>
{
    using ParentType = GridGeometryView<ActualGridGeometry>;
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
    { return this->gridGeometry_->numFaceDofs(); }
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
template<class GV, class T>
class StaggeredFVGridGeometry<GV, true, T>
: public BaseGridGeometry<GV, T>
{
    using ThisType = StaggeredFVGridGeometry<GV, true, T>;
    using ParentType = BaseGridGeometry<GV, T>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using Element = typename GV::template Codim<0>::Entity;

    using IntersectionMapper = typename T::IntersectionMapper;
    using GeometryHelper = typename T::GeometryHelper;
    using ConnectivityMap = typename T::template ConnectivityMap<ThisType>;

public:
    //! export the traits
    using Traits = typename T::PublicTraits;

    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::staggered;
    static constexpr int upwindSchemeOrder = T::upwindSchemeOrder;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;
    static constexpr bool cachingEnabled = true;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename T::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename T::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename T::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<T>;
    //! export the grid view type
    using GridView = GV;
    //! export the dof type indices
    using DofTypeIndices = typename T::DofTypeIndices;

    //! return a integral constant for cell center dofs
    static constexpr auto cellCenterIdx()
    { return typename DofTypeIndices::CellCenterIdx{}; }

    //! return a integral constant for face dofs
    static constexpr auto faceIdx()
    { return typename DofTypeIndices::FaceIdx{}; }

    //! The order of the stencil built
    static constexpr int upwindStencilOrder()
    {   return upwindSchemeOrder; }

    using CellCenterFVGridGeometryType = CellCenterFVGridGeometry<ThisType>;
    using FaceFVGridGeometryType = FaceFVGridGeometry<ThisType>;

    using FVGridGeometryTuple = std::tuple< CellCenterFVGridGeometry<ThisType>, FaceFVGridGeometry<ThisType> >;

    //! Constructor
    StaggeredFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
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
        scvfIndicesOfScv_.clear();
        intersectionMapper_.update();

        // determine size of containers
        std::size_t numScvs = this->gridView().size(0);
        std::size_t numScvf = 0;
        for (const auto& element : elements(this->gridView()))
            numScvf += element.subEntities(1);

        // reserve memory
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvf);
        scvfIndicesOfScv_.resize(numScvs);
        localToGlobalScvfIndices_.resize(numScvs);
        hasBoundaryScvf_.resize(numScvs, false);

        // Build the scvs and scv faces
        GridIndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            auto eIdx = this->elementMapper().index(element);

            // reserve memory for the localToGlobalScvfIdx map
            auto numLocalFaces = intersectionMapper_.numFaces(element);
            localToGlobalScvfIndices_[eIdx].resize(numLocalFaces);

            scvs_[eIdx] = SubControlVolume(element.geometry(), eIdx);

            // the element-wise index sets for finite volume geometry
            std::vector<GridIndexType> scvfsIndexSet;
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
                                        std::vector<GridIndexType>({eIdx, nIdx}),
                                        geometryHelper);
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;
                    scvfsIndexSet.push_back(scvfIdx++);
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    scvfs_.emplace_back(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<GridIndexType>({eIdx, this->gridView().size(0) + numBoundaryScvf_++}),
                                        geometryHelper);
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;
                    scvfsIndexSet.push_back(scvfIdx++);

                    hasBoundaryScvf_[eIdx] = true;
                }
            }

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
        }

        // build the connectivity map for an effecient assembly
        connectivityMap_.update(*this);
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        return scvfs_[scvfIdx];
    }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<GridIndexType>& scvfIndicesOfScv(GridIndexType scvIdx) const
    {
        return scvfIndicesOfScv_[scvIdx];
    }

    GridIndexType localToGlobalScvfIndex(GridIndexType eIdx, LocalIndexType localScvfIdx) const
    {
        return localToGlobalScvfIndices_[eIdx][localScvfIdx];
    }

    const SubControlVolumeFace& scvf(GridIndexType eIdx, LocalIndexType localScvfIdx) const
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

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf(GridIndexType eIdx) const
    { return hasBoundaryScvf_[eIdx]; }

private:

    // mappers
    ConnectivityMap connectivityMap_;
    IntersectionMapper intersectionMapper_;

    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<std::vector<GridIndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<GridIndexType>> localToGlobalScvfIndices_;
    GridIndexType numBoundaryScvf_;
    std::vector<bool> hasBoundaryScvf_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored.
 */
template<class GV, class T>
class StaggeredFVGridGeometry<GV, false, T>
: public BaseGridGeometry<GV, T>
{
    using ThisType = StaggeredFVGridGeometry<GV, false, T>;
    using ParentType = BaseGridGeometry<GV, T>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using Element = typename GV::template Codim<0>::Entity;

    using IntersectionMapper = typename T::IntersectionMapper;
    using ConnectivityMap = typename T::template ConnectivityMap<ThisType>;

public:
    //! export the traits
    using Traits = typename T::PublicTraits;

    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::staggered;
    static constexpr int upwindSchemeOrder = T::upwindSchemeOrder;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;
    static constexpr bool cachingEnabled = false;

    using GeometryHelper = typename T::GeometryHelper;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename T::template LocalView<ThisType, false>;
    //! export the type of sub control volume
    using SubControlVolume = typename T::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename T::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<T>;
    //! export the grid view type
    using GridView = GV;
    //! export the dof type indices
    using DofTypeIndices = typename T::DofTypeIndices;

    //! return a integral constant for cell center dofs
    static constexpr auto cellCenterIdx()
    { return typename DofTypeIndices::CellCenterIdx{}; }

    //! return a integral constant for face dofs
    static constexpr auto faceIdx()
    { return typename DofTypeIndices::FaceIdx{}; }

    //! The order of the stencil built
    static constexpr int upwindStencilOrder()
    {   return upwindSchemeOrder; }

    using CellCenterFVGridGeometryType = CellCenterFVGridGeometry<ThisType>;
    using FaceFVGridGeometryType = FaceFVGridGeometry<ThisType>;

    using FVGridGeometryTuple = std::tuple< CellCenterFVGridGeometry<ThisType>, FaceFVGridGeometry<ThisType> >;

    //! Constructor
    StaggeredFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
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
            std::vector<GridIndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);
            localToGlobalScvfIndices_[eIdx].resize(numLocalFaces);

            std::vector<GridIndexType> neighborVolVarIndexSet;
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

    const std::vector<GridIndexType>& scvfIndicesOfScv(GridIndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    GridIndexType localToGlobalScvfIndex(GridIndexType eIdx, LocalIndexType localScvfIdx) const
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
    const std::vector<GridIndexType>& neighborVolVarIndices(GridIndexType scvIdx) const
    { return neighborVolVarIndices_[scvIdx]; }

private:

    //! Information on the global number of geometries
    std::size_t numScvs_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;
    std::vector<std::vector<GridIndexType>> localToGlobalScvfIndices_;
    std::vector<std::vector<GridIndexType>> neighborVolVarIndices_;

    // mappers
    ConnectivityMap connectivityMap_;
    IntersectionMapper intersectionMapper_;

    //! vectors that store the global data
    std::vector<std::vector<GridIndexType>> scvfIndicesOfScv_;
};

} // end namespace

#endif
