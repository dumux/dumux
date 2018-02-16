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

    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::staggered;
    using LocalView = typename ActualFVGridGeometry::LocalView;

    static constexpr bool isCellCenter = false;
    static constexpr bool isFace = false;

    const auto& gridView() const
    { return fvGridGeometry_->gridView(); }

    const auto &connectivityMap() const // TODO return correct map
    { return fvGridGeometry_->connectivityMap(); }

    const auto& vertexMapper() const
    { return fvGridGeometry_->vertexMapper(); }

    const auto& elementMapper() const
    { return fvGridGeometry_->elementMapper(); }

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

    static constexpr bool isCellCenter = true;

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

    static constexpr bool isFace = true;

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
          bool enableFVGridGeometryCache,
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
    using Element = typename GV::template Codim<0>::Entity;

    enum {
        // Grid and world dimension
        dim = GV::dimension,
        dimWorld = GV::dimensionworld
    };

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

    using CellCenterFVGridGeometryType = CellCenterFVGridGeometry<ThisType>;
    using FaceFVGridGeometryType = FaceFVGridGeometry<ThisType>;

    using FVGridGeometryTuple = std::tuple< CellCenterFVGridGeometry<ThisType>, FaceFVGridGeometry<ThisType> >;

    //! Constructor
    StaggeredFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    , intersectionMapper_(gridView) {}

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

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return this->elementMap()[scv.elementIndex()]; }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return this->elementMap()[eIdx]; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update()
    {
        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
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
                                        geometryHelper
                                        );
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;
                    scvfsIndexSet.push_back(scvfIdx++);
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    scvfs_.emplace_back(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<IndexType>({eIdx, this->gridView().size(0) + numBoundaryScvf_++}),
                                        geometryHelper
                                        );
                    localToGlobalScvfIndices_[eIdx][localFaceIndex] = scvfIdx;
                    scvfsIndexSet.push_back(scvfIdx++);
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
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<IndexType>> localToGlobalScvfIndices_;
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
{
    // TODO: implement without caching
};

} // end namespace

#endif
