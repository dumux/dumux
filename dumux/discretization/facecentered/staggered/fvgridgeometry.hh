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
 * \ingroup FaceCenteredStaggeredFVElementGeometry
 * \copydoc Dumux::FaceCenteredStaggeredFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_GRID_GEOMETRY
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_GRID_GEOMETRY

#include <memory>

#include <dumux/common/defaultmappertraits.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/math.hh>

#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/discretization/facecentered/staggered/subcontrolvolume.hh>
#include <dumux/discretization/facecentered/staggered/subcontrolvolumeface.hh>
#include <dumux/discretization/facecentered/staggered/fvelementgeometry.hh>
#include <dumux/discretization/facecentered/staggered/geometryhelper.hh>
#include <dumux/discretization/facecentered/staggered/connectivitymap.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief The default traits for the face-center staggered finite volume grid geometry
 *        Defines the scv and scvf types and the mapper types
 * \tparam GridView the grid view type
 */
template<class GridView>
struct FaceCenteredStaggeredDefaultGridGeometryTraits : public DefaultMapperTraits<GridView>
{
    using SubControlVolume = FaceCenteredStaggeredSubControlVolume<GridView>;
    using SubControlVolumeFace = FaceCenteredStaggeredSubControlVolumeFace<GridView>;
    using IntersectionMapper = ConformingGridIntersectionMapper<GridView>;
    using GeometryHelper = FaceCenteredStaggeredGeometryHelper<GridView, typename GridView::Grid>;

    template<class GridGeometry>
    using ConnectivityMap = FaceCenteredStaggeredConnectivityMap<GridGeometry>;

    template<class GridGeometry, bool enableCache>
    using LocalView = FaceCenteredStaggeredFVElementGeometry<GridGeometry, enableCache>;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for face-centered staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class GridView,
         bool cachingEnabled = false,
         class Traits = FaceCenteredStaggeredDefaultGridGeometryTraits<GridView>>
class FaceCenteredStaggeredFVGridGeometry;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored.
 */
template<class GV, class Traits>
class FaceCenteredStaggeredFVGridGeometry<GV, true, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = FaceCenteredStaggeredFVGridGeometry<GV, true, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using SmallLocalIndexType = typename IndexTraits<GV>::SmallLocalIndex;
    using Element = typename GV::template Codim<0>::Entity;

    using IntersectionMapper = typename Traits::IntersectionMapper;
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;

    using Scalar = typename GV::ctype;

    static constexpr auto dim = GV::Grid::dimension;
    static constexpr auto numFacesPerElement = dim * 2;
    static constexpr auto numScvsPerElement = numFacesPerElement;
    static constexpr auto numLateralScvfsPerScv = 2 * (dim - 1);
    static constexpr auto numLateralScvfsPerElement = numFacesPerElement*numLateralScvfsPerScv;

    static constexpr auto maxNumScvfsPerElement = numLateralScvfsPerElement  // number of lateral faces
                                                + numFacesPerElement  // number of central frontal faces
                                                + numFacesPerElement; // number of potential frontal faces on boundary

    // static_assert(maxNumScvfsPerElement == 16); // TODO remove

    using ScvfCornerStorage = typename Traits::SubControlVolumeFace::Traits::CornerStorage;
    using ScvCornerStorage = typename Traits::SubControlVolume::Traits::CornerStorage;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::FCStaggered;
    static constexpr DiscretizationMethod discMethod{};

    static constexpr bool cachingEnabled = true;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the grid view type
    using GridView = GV;
    //! export the geometry helper type
    using GeometryHelper = typename Traits::GeometryHelper;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;

    //! Constructor
    FaceCenteredStaggeredFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
    : ParentType(gridView)
    , intersectionMapper_(gridView)
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethods::fcstaggered>::isValid(gridView))
            DUNE_THROW(Dune::InvalidStateException, "The staggered discretization method needs at least an overlap of 1 for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");

        update_();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return scvs_.size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return scvfs_.size(); }

    //! The total number of boundary sub control volumes
    std::size_t numBoundaryScv() const
    { return numBoundaryScv_; }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! The total number of intersections
    std::size_t numIntersections() const
    { return intersectionMapper_.numIntersections(); }

    //! the total number of dofs
    std::size_t numDofs() const
    { return this->gridView().size(1); }

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

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

    //! Get the global sub control volume indices of an element
    const std::array<GridIndexType, numScvsPerElement>& scvIndicesOfElement(GridIndexType eIdx) const
    { return scvIndicesOfElement_[eIdx]; }

    //! Get the global sub control volume face indices of an element
    const std::vector<GridIndexType>& scvfIndicesOfElement(GridIndexType eIdx) const
    { return scvfIndicesOfElement_[eIdx]; }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const ConnectivityMap& connectivityMap() const
    { return connectivityMap_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf(GridIndexType eIdx) const
    { return hasBoundaryScvf_[eIdx]; }

    //! Return a reference to the intersection mapper
    const IntersectionMapper& intersectionMapper() const
    { return intersectionMapper_; }

    SmallLocalIndexType firstLocalScvfIdxOfScv(const GridIndexType eIdx, const SmallLocalIndexType localScvIdx) const
    { return scvfOfScvInfo_[eIdx][localScvIdx]; }

    //! If a d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicFaceMap_.count(dofIdx); }

    //! The index of the d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicFaceMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries // TODO rename to periodic dof map in fvassembler
    const std::unordered_map<GridIndexType, GridIndexType>& periodicVertexMap() const
    { return periodicFaceMap_; }

    //! get the global index of the orthogonal face sharing a common entity
    GridIndexType lateralOrthogonalScvf(const SubControlVolumeFace& scvf) const
    { return lateralOrthogonalScvf_[scvf.index()]; }

private:

    void update_()
    {
        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvIndicesOfElement_.clear();
        scvfIndicesOfElement_.clear();
        lateralOrthogonalScvf_.clear();
        intersectionMapper_.update(this->gridView());

        // determine size of containers
        const auto numElements = this->gridView().size(0);
        scvIndicesOfElement_.resize(numElements);
        scvfIndicesOfElement_.resize(numElements);
        hasBoundaryScvf_.resize(numElements, false);
        scvfOfScvInfo_.resize(numElements);
        // on frontal + maybe boundary per face and  2*(dim-1) laterals per frontal
        lateralOrthogonalScvf_.resize(numElements*(2*dim*(2 + 2*(dim-1))));

        outSideBoundaryVolVarIdx_ = 0;
        numBoundaryScv_ = 0;
        numBoundaryScvf_ = 0;

        GeometryHelper geometryHelper(this->gridView());

        // get the global scv indices first
        GridIndexType scvIdx = 0;
        GridIndexType scvfIdx = 0;
        for (const auto& element : elements(this->gridView()))
        {
            assert(numScvsPerElement == element.subEntities(1));

            // the element-wise index sets for finite volume geometry
            std::array<GridIndexType, numScvsPerElement> scvsIndexSet;

            std::vector<GridIndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(1/*frontal in element*/ + 1 /*frontal on boundary*/ + numLateralScvfsPerScv);

            // keep track of frontal boundary scvfs
            std::size_t numFrontalBoundaryScvfs = 0;

            // a temporary map to store pairs of common entities
            std::unordered_map<GridIndexType, Dune::ReservedVector<GridIndexType, 2>> commonEntityIdxToScvfsMap;

            SmallLocalIndexType localScvIdx = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                scvsIndexSet[localScvIdx] = scvIdx++; // one scv per element face

                // the frontal sub control volume face at the element center
                scvfsIndexSet.push_back(scvfIdx++);

                if constexpr(dim > 1)
                {
                    // the lateral sub control volume faces
                    for (const auto lateralFacetIndex : geometryHelper.localLaterFaceIndices(localScvIdx))
                    {
                        const auto& lateralIntersection = geometryHelper.intersection(lateralFacetIndex, element);
                        if (onProcessorBoundary_(lateralIntersection))
                            continue;

                        scvfsIndexSet.push_back(scvfIdx);

                        const auto commonEntityIdx = geometryHelper.globalCommonEntityIndex(
                            element, localScvIdx, lateralFacetIndex
                        );
                        commonEntityIdxToScvfsMap[commonEntityIdx].push_back(scvfIdx);
                        ++scvfIdx;
                    }
                }

                // handle physical domain boundary
                if (onDomainBoundary_(intersection))
                {
                    ++numBoundaryScv_; // frontal face
                    numBoundaryScv_ += numLateralScvfsPerScv; // boundary scvs for lateral faces
                    ++numFrontalBoundaryScvfs;
                }

                ++localScvIdx;
            }

            // add global indices of frontal boundary scvfs last
            for (std::size_t i = 0; i < numFrontalBoundaryScvfs; ++i)
                scvfsIndexSet.push_back(scvfIdx++);

            // Save the scv indices belonging to this element to build up fv element geometries fast
            const auto eIdx = this->elementMapper().index(element);
            scvIndicesOfElement_[eIdx] = std::move(scvsIndexSet);
            scvfIndicesOfElement_[eIdx] = std::move(scvfsIndexSet);

            // create the bi-directional map
            for (const auto& [_, scvfs] : commonEntityIdxToScvfsMap)
            {
                // TODO: this maybe be less than 2 if there is a processor boundary
                // are we sure that the lateral orthogonal scvf is then never needed?
                if (scvfs.size() == 2)
                {
                    lateralOrthogonalScvf_[scvfs[0]] = scvfs[1];
                    lateralOrthogonalScvf_[scvfs[1]] = scvfs[0];
                }
            }
        }

         // reserve memory
        const auto numInteriorScvs = scvIdx;
        scvs_.reserve(numInteriorScvs);
        scvfs_.reserve((numScvsPerElement/*one interior frontal scvf per scv*/
                        +numLateralScvfsPerElement)*numElements + numBoundaryScv_);

        // Build the scvs and scv faces
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);
            const auto& globalScvIndices = scvIndicesOfElement_[eIdx];
            const auto& globalScvfIndices = scvfIndicesOfElement_[eIdx];
            scvfOfScvInfo_[eIdx].reserve(globalScvfIndices.size());

            SmallLocalIndexType localScvIdx = 0; // also corresponds to local element face index (one scv per face)
            SmallLocalIndexType localScvfIdx = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                const auto dofIndex = intersectionMapper().globalIntersectionIndex(element, localScvIdx);
                const auto localOppositeScvIdx = geometryHelper.localOppositeIdx(localScvIdx);
                const auto& intersectionGeometry = intersection.geometry();
                const auto& elementGeometry = element.geometry();

                // store the local index of the first scvf corresponding to the scv on the current intersection
                scvfOfScvInfo_[eIdx].push_back(localScvfIdx);

                // handle periodic boundaries
                if (onPeriodicBoundary_(intersection))
                {
                    this->setPeriodic();

                    const auto& otherElement = intersection.outside();

                    SmallLocalIndexType otherIntersectionLocalIdx = 0;
                    bool periodicFaceFound = false;

                    for (const auto& otherIntersection : intersections(this->gridView(), otherElement))
                    {
                        if (periodicFaceFound)
                            continue;

                        if (Dune::FloatCmp::eq(intersection.centerUnitOuterNormal()*otherIntersection.centerUnitOuterNormal(), -1.0, 1e-7))
                        {
                            const auto periodicDofIdx = intersectionMapper().globalIntersectionIndex(otherElement, otherIntersectionLocalIdx);
                            periodicFaceMap_[dofIndex] = periodicDofIdx;
                            periodicFaceFound = true;
                        }

                        ++otherIntersectionLocalIdx;
                    }
                }

                // the sub control volume
                scvs_.emplace_back(elementGeometry,
                                   intersectionGeometry,
                                   globalScvIndices[localScvIdx],
                                   localScvIdx,
                                   dofIndex,
                                   Dumux::normalAxis(intersection.centerUnitOuterNormal()),
                                   this->elementMapper().index(element),
                                   onDomainBoundary_(intersection));

                // the frontal sub control volume face at the element center
                scvfs_.emplace_back(elementGeometry,
                                    intersectionGeometry,
                                    std::array{globalScvIndices[localScvIdx], globalScvIndices[localOppositeScvIdx]},
                                    localScvfIdx,
                                    globalScvfIndices[localScvfIdx],
                                    intersection.centerUnitOuterNormal(),
                                    SubControlVolumeFace::FaceType::frontal,
                                    false);
                ++localScvfIdx;

                // the lateral sub control volume faces
                const auto lateralFaceIndices = geometryHelper.localLaterFaceIndices(localScvIdx);
                for (const auto lateralFacetIndex : lateralFaceIndices)
                {
                    const auto& lateralIntersection = geometryHelper.intersection(lateralFacetIndex, element);
                    if (onProcessorBoundary_(lateralIntersection))
                        continue;

                    const auto& lateralFacet =  geometryHelper.facet(lateralFacetIndex, element);
                    const auto& lateralFacetGeometry = lateralFacet.geometry();

                    if (lateralIntersection.neighbor()) // TODO: periodic?
                        geometryHelper.update(element, lateralIntersection.outside());

                    // helper lambda to get the lateral scvf's global inside and outside scv indices
                    const auto globalScvIndicesForLateralFace = [&]
                    {
                        const auto globalOutsideScvIdx = [&]
                        {
                            if (lateralIntersection.neighbor())
                            {
                                const auto parallelElemIdx = this->elementMapper().index(lateralIntersection.outside());
                                const auto localOutsideScvIdx = geometryHelper.localFaceIndexInOtherElement(localScvIdx);
                                const auto& globalScvIndicesOfNeighborElement = scvIndicesOfElement_[parallelElemIdx];
                                return globalScvIndicesOfNeighborElement[localOutsideScvIdx];
                            }
                            else
                                return numInteriorScvs + outSideBoundaryVolVarIdx_++;
                        }();

                        return std::array{globalScvIndices[localScvIdx], globalOutsideScvIdx};
                    }();

                    scvfs_.emplace_back(elementGeometry,
                                        intersectionGeometry,
                                        lateralFacetGeometry,
                                        globalScvIndicesForLateralFace, // TODO higher order
                                        localScvfIdx,
                                        globalScvfIndices[localScvfIdx],
                                        lateralIntersection.centerUnitOuterNormal(),
                                        SubControlVolumeFace::FaceType::lateral,
                                        onDomainBoundary_(lateralIntersection));
                    ++localScvfIdx;

                    if (onDomainBoundary_(lateralIntersection))
                    {
                        ++numBoundaryScvf_;
                        hasBoundaryScvf_[eIdx] = true;
                    }
                } // end loop over lateral facets

                ++localScvIdx;
            } // end first loop over intersections

            // do a second loop over all intersections to add frontal boundary faces
            localScvIdx = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                // the frontal sub control volume face at a domain boundary (coincides with element face)
                if (onDomainBoundary_(intersection))
                {
                    ++numBoundaryScvf_;

                    // the frontal sub control volume face at the boundary
                    scvfs_.emplace_back(element.geometry(),
                                        intersection.geometry(),
                                        std::array{globalScvIndices[localScvIdx], globalScvIndices[localScvIdx]}, // TODO outside boundary, periodic, parallel?
                                        localScvfIdx,
                                        globalScvfIndices[localScvfIdx],
                                        intersection.centerUnitOuterNormal(),
                                        SubControlVolumeFace::FaceType::frontal,
                                        true);
                    ++localScvfIdx;
                    hasBoundaryScvf_[eIdx] = true;
                }

                ++localScvIdx;
            }
        }

        connectivityMap_.update(*this);

        // allow to free the non-needed memory
        lateralOrthogonalScvf_.resize(scvfs_.size());
    }

    bool onDomainBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && intersection.boundary();
    }

    bool onProcessorBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && !intersection.boundary();
    }

    bool onPeriodicBoundary_(const typename GridView::Intersection& intersection) const
    {
        return intersection.boundary() && intersection.neighbor();
    }

    // mappers
    ConnectivityMap connectivityMap_;
    IntersectionMapper intersectionMapper_;

    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    GridIndexType numBoundaryScv_;
    GridIndexType numBoundaryScvf_;
    GridIndexType outSideBoundaryVolVarIdx_;
    std::vector<bool> hasBoundaryScvf_;
    std::vector<std::vector<SmallLocalIndexType>> scvfOfScvInfo_;

    std::vector<GridIndexType> lateralOrthogonalScvf_;
    std::vector<std::array<GridIndexType, numScvsPerElement>> scvIndicesOfElement_;
    std::vector<std::vector<GridIndexType>> scvfIndicesOfElement_;

    // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicFaceMap_;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for face-centered staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored.
 */
template<class GV, class Traits>
class FaceCenteredStaggeredFVGridGeometry<GV, false, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = FaceCenteredStaggeredFVGridGeometry<GV, false, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using SmallLocalIndexType = typename IndexTraits<GV>::SmallLocalIndex;
    using Element = typename GV::template Codim<0>::Entity;

    using IntersectionMapper = typename Traits::IntersectionMapper;
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;

    static constexpr auto dim = GV::Grid::dimension;
    static constexpr auto numFacesPerElement = dim * 2;
    static constexpr auto numScvsPerElement = numFacesPerElement;
    static constexpr auto numLateralScvfsPerScv = 2 * (dim - 1);
    static constexpr auto numLateralScvfsPerElement = numFacesPerElement*numLateralScvfsPerScv;

    static constexpr auto maxNumScvfsPerElement = numLateralScvfsPerElement  // number of lateral faces
                                                + numFacesPerElement  // number of central frontal faces
                                                + numFacesPerElement; // number of potential frontal faces on boundary

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::FCStaggered;
    static constexpr DiscretizationMethod discMethod{};

    static constexpr bool cachingEnabled = false;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, false>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the grid view type
    using GridView = GV;
    //! export the geometry helper type
    using GeometryHelper = typename Traits::GeometryHelper;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;

    //! Constructor
    FaceCenteredStaggeredFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
    : ParentType(gridView)
    , intersectionMapper_(gridView)
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethods::fcstaggered>::isValid(gridView))
            DUNE_THROW(Dune::InvalidStateException, "The staggered discretization method needs at least an overlap of 1 for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");

        update_();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return numScvs_; }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return numScvf_; }

    //! The total number of boundary sub control volumes
    std::size_t numBoundaryScv() const
    { return numBoundaryScv_; }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! The total number of intersections
    std::size_t numIntersections() const
    { return intersectionMapper_.numIntersections(); }

    //! the total number of dofs
    std::size_t numDofs() const
    { return this->gridView().size(1); }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const ConnectivityMap& connectivityMap() const
    { return connectivityMap_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf(GridIndexType eIdx) const
    { return hasBoundaryScvf_[eIdx]; }

    //! Return a reference to the intersection mapper
    const IntersectionMapper& intersectionMapper() const
    { return intersectionMapper_; }

    //! Get the global sub control volume face indices of an element
    const std::vector<GridIndexType>& scvfIndicesOfElement(GridIndexType eIdx) const
    { return scvfIndicesOfElement_[eIdx]; }

    //! Get the global sub control volume face indices of an element
    GridIndexType outsideVolVarIndex(GridIndexType scvfIdx) const
    { return outsideVolVarIndices_.at(scvfIdx); }

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

    //! If a d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicFaceMap_.count(dofIdx); }

    //! The index of the d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicFaceMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries // TODO rename to periodic dof map in fvassembler
    const std::unordered_map<GridIndexType, GridIndexType>& periodicVertexMap() const
    { return periodicFaceMap_; }

    //! get the global index of the orthogonal face sharing a common entity
    GridIndexType lateralOrthogonalScvf(const SubControlVolumeFace& scvf) const
    { return lateralOrthogonalScvf_[scvf.index()]; }

private:

    void update_()
    {
        intersectionMapper_.update(this->gridView());

        // clear local data
        numScvs_ = 0;
        numScvf_ = 0;
        numBoundaryScv_ = 0;
        numBoundaryScvf_ = 0;
        hasBoundaryScvf_.clear();
        scvfIndicesOfElement_.clear();
        lateralOrthogonalScvf_.clear();
        outsideVolVarIndices_.clear();

        // determine size of containers
        const auto numElements = this->gridView().size(0);
        scvfIndicesOfElement_.resize(numElements);
        hasBoundaryScvf_.resize(numElements, false);
        // on frontal + maybe boundary per face and  2*(dim-1) laterals per frontal
        lateralOrthogonalScvf_.resize(numElements*(2*dim*(2 + 2*(dim-1))));

        GeometryHelper geometryHelper(this->gridView());

        // get the global scv indices first
        GridIndexType scvfIdx = 0;

        GridIndexType neighborVolVarIdx = numElements*numScvsPerElement;

        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);
            assert(numScvsPerElement == element.subEntities(1));

            // the element-wise index sets for finite volume geometry
            std::vector<GridIndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(1/*frontal in element*/ + 1 /*frontal on boundary*/ + numLateralScvfsPerScv);

            // keep track of frontal boundary scvfs
            std::size_t numFrontalBoundaryScvfs = 0;

            // a temporary map to store pairs of common entities
            std::unordered_map<GridIndexType, Dune::ReservedVector<GridIndexType, 2>> commonEntityIdxToScvfsMap;

            SmallLocalIndexType localScvIdx = 0;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                ++numScvs_;

                // the frontal sub control volume face at the element center
                scvfsIndexSet.push_back(scvfIdx++);

                if constexpr(dim > 1)
                {
                    // the lateral sub control volume faces
                    for (const auto lateralFacetIndex : geometryHelper.localLaterFaceIndices(localScvIdx))
                    {
                        const auto& lateralIntersection = geometryHelper.intersection(lateralFacetIndex, element);
                        if (onProcessorBoundary_(lateralIntersection))
                            continue;

                        if (onDomainBoundary_(lateralIntersection))
                            outsideVolVarIndices_[scvfIdx] = neighborVolVarIdx++;

                        scvfsIndexSet.push_back(scvfIdx);

                        const auto commonEntityIdx = geometryHelper.globalCommonEntityIndex(
                            element, localScvIdx, lateralFacetIndex
                        );
                        commonEntityIdxToScvfsMap[commonEntityIdx].push_back(scvfIdx);

                        ++scvfIdx;
                    }
                }

                // handle physical domain boundary
                if (onDomainBoundary_(intersection))
                {
                    ++numBoundaryScv_; // frontal face
                    numBoundaryScv_ += numLateralScvfsPerScv; // boundary scvs for lateral faces
                    ++numFrontalBoundaryScvfs;
                    ++numBoundaryScvf_;
                    hasBoundaryScvf_[eIdx] = true;
                }

                // handle periodic boundaries
                if (onPeriodicBoundary_(intersection))
                {
                    this->setPeriodic();

                    const auto& otherElement = intersection.outside();

                    SmallLocalIndexType otherIntersectionLocalIdx = 0;
                    bool periodicFaceFound = false;

                    for (const auto& otherIntersection : intersections(this->gridView(), otherElement))
                    {
                        if (periodicFaceFound)
                            continue;

                        if (Dune::FloatCmp::eq(intersection.centerUnitOuterNormal()*otherIntersection.centerUnitOuterNormal(), -1.0, 1e-7))
                        {
                            const auto periodicDofIdx = intersectionMapper().globalIntersectionIndex(otherElement, otherIntersectionLocalIdx);
                            const auto dofIndex = intersectionMapper().globalIntersectionIndex(element, localScvIdx);
                            periodicFaceMap_[dofIndex] = periodicDofIdx;
                            periodicFaceFound = true;
                        }

                        ++otherIntersectionLocalIdx;
                    }
                }

                // the lateral sub control volume faces
                const auto lateralFaceIndices = geometryHelper.localLaterFaceIndices(localScvIdx);
                for (const auto lateralFacetIndex : lateralFaceIndices)
                {
                    const auto& lateralIntersection = geometryHelper.intersection(lateralFacetIndex, element);
                    if (onProcessorBoundary_(lateralIntersection))
                        continue;

                    if (onDomainBoundary_(lateralIntersection))
                    {
                        ++numBoundaryScvf_;
                        hasBoundaryScvf_[eIdx] = true;
                    }
                } // end loop over lateral facets

                ++localScvIdx;
            }

            // add global indices of frontal boundary scvfs last
            for (std::size_t i = 0; i < numFrontalBoundaryScvfs; ++i)
                scvfsIndexSet.push_back(scvfIdx++);

            // Save the scv indices belonging to this element to build up fv element geometries fast
            scvfIndicesOfElement_[eIdx] = std::move(scvfsIndexSet);

            // create the bi-directional map
            for (const auto& [_, scvfs] : commonEntityIdxToScvfsMap)
            {
                // TODO: this maybe be less than 2 if there is a processor boundary
                // are we sure that the lateral orthogonal scvf is then never needed?
                if (scvfs.size() == 2)
                {
                    lateralOrthogonalScvf_[scvfs[0]] = scvfs[1];
                    lateralOrthogonalScvf_[scvfs[1]] = scvfs[0];
                }
            }
        }

        // set number of subcontrolvolume faces
        numScvf_ = scvfIdx;

        connectivityMap_.update(*this);
    }

    bool onDomainBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && intersection.boundary();
    }

    bool onProcessorBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && !intersection.boundary();
    }

    bool onPeriodicBoundary_(const typename GridView::Intersection& intersection) const
    {
        return intersection.boundary() && intersection.neighbor();
    }

    // mappers
    ConnectivityMap connectivityMap_;
    IntersectionMapper intersectionMapper_;

    //! Information on the global number of geometries
    std::size_t numScvs_;
    std::size_t numScvf_;
    std::size_t numBoundaryScv_;
    std::size_t numBoundaryScvf_;
    std::vector<bool> hasBoundaryScvf_;

    std::vector<GridIndexType> lateralOrthogonalScvf_;
    std::vector<std::vector<GridIndexType>> scvfIndicesOfElement_;

    // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicFaceMap_;
    std::unordered_map<GridIndexType, GridIndexType> outsideVolVarIndices_;
};

} // end namespace Dumux

#endif
