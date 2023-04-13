// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_GRID_GEOMETRY
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_GRID_GEOMETRY

#include <memory>

#include <dune/common/rangeutilities.hh>
#include <dune/grid/common/scsgmapper.hh>

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
#include <dumux/discretization/facecentered/staggered/localintersectionindexmapper.hh>

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
    using LocalIntersectionMapper = FaceCenteredStaggeredLocalIntersectionIndexMapper<GridView>;
    using GeometryHelper = FaceCenteredStaggeredGeometryHelper<GridView>;

    template<class GridGeometry>
    using ConnectivityMap = FaceCenteredStaggeredConnectivityMap<GridGeometry>;

    template<class GridGeometry, bool enableCache>
    using LocalView = FaceCenteredStaggeredFVElementGeometry<GridGeometry, enableCache>;

    struct StaticInfo
    {
        static constexpr auto dim = GridView::Grid::dimension;
        static constexpr auto numFacesPerElement = dim * 2;
        static constexpr auto numScvsPerElement = numFacesPerElement;
        static constexpr auto numLateralScvfsPerScv = 2 * (dim - 1);
        static constexpr auto numLateralScvfsPerElement = numFacesPerElement*numLateralScvfsPerScv;
        static constexpr auto minNumScvfsPerElement = numLateralScvfsPerElement  // number of lateral faces
                                                    + numFacesPerElement;  // number of central frontal faces
        static constexpr auto maxNumScvfsPerElement = minNumScvfsPerElement
                                                    + numFacesPerElement; // number of potential frontal faces on boundary
    };
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
 * \ingroup FaceCenteredStaggeredDiscretization
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

    static constexpr auto dim = Traits::StaticInfo::dim;
    static constexpr auto numScvsPerElement = Traits::StaticInfo::numScvsPerElement;
    static constexpr auto numLateralScvfsPerScv = Traits::StaticInfo::numLateralScvfsPerScv;
    static constexpr auto numLateralScvfsPerElement = Traits::StaticInfo::numLateralScvfsPerElement;
    static constexpr auto minNumScvfsPerElement = Traits::StaticInfo::minNumScvfsPerElement;
    static constexpr auto maxNumScvfsPerElement = Traits::StaticInfo::maxNumScvfsPerElement;

    using ScvfCornerStorage = typename Traits::SubControlVolumeFace::Traits::CornerStorage;
    using ScvCornerStorage = typename Traits::SubControlVolume::Traits::CornerStorage;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::FCStaggered;
    static constexpr DiscretizationMethod discMethod{};

    static constexpr bool cachingEnabled = true;

    //! export basic grid geometry type for the alternative constructor
    using BasicGridGeometry = BasicGridGeometry_t<GV, Traits>;
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
    //! export the local intersection mapper
    using LocalIntersectionMapper = typename Traits::LocalIntersectionMapper;
    //! export static information
    using StaticInformation = typename Traits::StaticInfo;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;

    //! Constructor with basic grid geometry used to share state with another grid geometry on the same grid view
    FaceCenteredStaggeredFVGridGeometry(std::shared_ptr<BasicGridGeometry> gg, const std::string& paramGroup = "")
    : ParentType(std::move(gg))
    , intersectionMapper_(this->gridView())
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod>::isValid(this->gridView()))
            DUNE_THROW(Dune::InvalidStateException, "The staggered discretization method needs at least an overlap of 1 for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");

        update_();
    }

    //! Constructor from gridView
    FaceCenteredStaggeredFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
    : FaceCenteredStaggeredFVGridGeometry(std::make_shared<BasicGridGeometry>(gridView), paramGroup)
    {}

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

    //! Iterator range for sub control volumes. Iterates over
    //! all scvs of the element-local fvGeometry.
    auto scvs(const LocalView& fvGeometry) const
    {
        auto begin = scvs_.cbegin() + numScvsPerElement*fvGeometry.elementIndex();
        const auto end = begin + numScvsPerElement;
        return Dune::IteratorRange<std::decay_t<decltype(begin)>>(begin, end);
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

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

    //! If a d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { return periodicFaceMap_.count(dofIdx); }

    //! The index of the d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { return periodicFaceMap_.at(dofIdx); }

    //! Returns the map between dofs across periodic boundaries // TODO rename to periodic dof map in fvassembler
    const std::unordered_map<GridIndexType, GridIndexType>& periodicVertexMap() const
    { return periodicFaceMap_; }

private:

    void update_()
    {
        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvfIndicesOfElement_.clear();
        intersectionMapper_.update(this->gridView());

        // determine size of containers
        const auto numElements = this->gridView().size(0);
        scvfIndicesOfElement_.resize(numElements);
        hasBoundaryScvf_.resize(numElements, false);

        outSideBoundaryVolVarIdx_ = 0;
        numBoundaryScv_ = 0;
        numBoundaryScvf_ = 0;

        GeometryHelper geometryHelper(this->gridView());

        // get the global scvf indices first
        GridIndexType numScvfs = 0;
        for (const auto& element : elements(this->gridView()))
        {
            assert(numScvsPerElement == element.subEntities(1));

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                // frontal scvf in element center
                ++numScvfs;

                // lateral scvfs
                numScvfs += numLateralScvfsPerScv;

                // handle physical domain boundary
                if (onDomainBoundary_(intersection))
                {
                    ++numBoundaryScv_; // frontal face
                    numBoundaryScv_ += numLateralScvfsPerScv; // boundary scvs for lateral faces

                    // frontal scvf at boundary
                    ++numScvfs;
                }
            }
        }

         // allocate memory
        const auto numScvs = numElements*numScvsPerElement;
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvfs);

        // Build the scvs and scv faces
        std::size_t globalScvfIdx = 0;
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);
            auto& globalScvfIndices = scvfIndicesOfElement_[eIdx];
            globalScvfIndices.resize(minNumScvfsPerElement);
            globalScvfIndices.reserve(maxNumScvfsPerElement);

            auto getGlobalScvIdx = [&](const auto elementIdx, const auto localScvIdx)
            { return numScvsPerElement*elementIdx + localScvIdx; };

            LocalIntersectionMapper localIsMapper;
            localIsMapper.update(this->gridView(), element);

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                const auto& intersectionUnitOuterNormal = intersection.centerUnitOuterNormal();
                const auto localScvIdx = localIsMapper.realToRefIdx(intersection.indexInInside());
                auto localScvfIdx = localScvIdx*(1 + numLateralScvfsPerScv);

                const auto globalScvIdx = getGlobalScvIdx(eIdx, localScvIdx);
                const auto dofIndex = intersectionMapper().globalIntersectionIndex(element, intersection.indexInInside());
                const auto localOppositeScvIdx = geometryHelper.localOppositeIdx(localScvIdx);
                const auto& intersectionGeometry = intersection.geometry();
                const auto& elementGeometry = element.geometry();

                assert(localIsMapper.refToRealIdx(localScvIdx) == intersection.indexInInside());

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

                        if (Dune::FloatCmp::eq(intersectionUnitOuterNormal*otherIntersection.centerUnitOuterNormal(), -1.0, 1e-7))
                        {
                            const auto periodicDofIdx = intersectionMapper().globalIntersectionIndex(otherElement, otherIntersectionLocalIdx);
                            periodicFaceMap_[dofIndex] = periodicDofIdx;
                            periodicFaceFound = true;
                        }

                        ++otherIntersectionLocalIdx;
                    }
                }

                // the sub control volume
                scvs_[globalScvIdx] = SubControlVolume(
                    elementGeometry,
                    intersectionGeometry,
                    globalScvIdx,
                    localScvIdx,
                    dofIndex,
                    Dumux::normalAxis(intersectionUnitOuterNormal),
                    this->elementMapper().index(element),
                    onDomainBoundary_(intersection)
                 );

                // the frontal sub control volume face at the element center
                scvfs_.emplace_back(elementGeometry,
                    intersectionGeometry,
                    std::array{globalScvIdx, getGlobalScvIdx(eIdx, localOppositeScvIdx)},
                    localScvfIdx,
                    globalScvfIdx,
                    intersectionUnitOuterNormal,
                    SubControlVolumeFace::FaceType::frontal,
                    SubControlVolumeFace::BoundaryType::interior
                );

                globalScvfIndices[localScvfIdx] = globalScvfIdx++;
                ++localScvfIdx;

                // the lateral sub control volume faces
                for (const auto lateralFacetIndex : Dune::transformedRangeView(geometryHelper.localLaterFaceIndices(localScvIdx),
                                                                                [&](auto&& idx) { return localIsMapper.refToRealIdx(idx) ;})
                    )
                {
                    const auto& lateralIntersection = geometryHelper.intersection(lateralFacetIndex, element);

                    // helper lambda to get the lateral scvf's global inside and outside scv indices
                    const auto globalScvIndicesForLateralFace = [&]
                    {
                        const auto globalOutsideScvIdx = [&]
                        {
                            if (lateralIntersection.neighbor())
                            {
                                const auto parallelElemIdx = this->elementMapper().index(lateralIntersection.outside());
                                return getGlobalScvIdx(parallelElemIdx, localScvIdx);
                            }
                            else if (onDomainBoundary_(lateralIntersection))
                                return numScvs + outSideBoundaryVolVarIdx_++;
                            else
                                return globalScvIdx; // fallback for parallel, won't be used anyway
                        }();

                        return std::array{globalScvIdx, globalOutsideScvIdx};
                    }();

                    const auto boundaryType = [&]
                    {
                        if (onProcessorBoundary_(lateralIntersection))
                            return SubControlVolumeFace::BoundaryType::processorBoundary;
                        else if (onDomainBoundary_(lateralIntersection))
                            return SubControlVolumeFace::BoundaryType::physicalBoundary;
                        else
                            return SubControlVolumeFace::BoundaryType::interior;
                    }();

                    scvfs_.emplace_back(
                        elementGeometry,
                        intersectionGeometry,
                        geometryHelper.facet(lateralFacetIndex, element).geometry(),
                        globalScvIndicesForLateralFace, // TODO higher order
                        localScvfIdx,
                        globalScvfIdx,
                        lateralIntersection.centerUnitOuterNormal(),
                        SubControlVolumeFace::FaceType::lateral,
                        boundaryType
                    );

                    globalScvfIndices[localScvfIdx] = globalScvfIdx++;
                    ++localScvfIdx;

                    if (onDomainBoundary_(lateralIntersection))
                    {
                        ++numBoundaryScvf_;
                        hasBoundaryScvf_[eIdx] = true;
                    }
                } // end loop over lateral facets

            } // end first loop over intersections

            // do a second loop over all intersections to add frontal boundary faces
            int localScvfIdx = minNumScvfsPerElement;
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                // the frontal sub control volume face at a domain boundary (coincides with element face)
                if (onDomainBoundary_(intersection))
                {
                    const auto localScvIdx = localIsMapper.realToRefIdx(intersection.indexInInside());
                    const auto globalScvIdx = getGlobalScvIdx(eIdx, localScvIdx);
                    ++numBoundaryScvf_;

                    // the frontal sub control volume face at the boundary
                    scvfs_.emplace_back(
                        element.geometry(),
                        intersection.geometry(),
                        std::array{globalScvIdx, globalScvIdx}, // TODO outside boundary, periodic, parallel?
                        localScvfIdx,
                        globalScvfIdx,
                        intersection.centerUnitOuterNormal(),
                        SubControlVolumeFace::FaceType::frontal,
                        SubControlVolumeFace::BoundaryType::physicalBoundary
                    );

                    globalScvfIndices.push_back(globalScvfIdx);
                    ++globalScvfIdx;
                    ++localScvfIdx;
                    hasBoundaryScvf_[eIdx] = true;
                }
            }
        }

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

    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    GridIndexType numBoundaryScv_;
    GridIndexType numBoundaryScvf_;
    GridIndexType outSideBoundaryVolVarIdx_;
    std::vector<bool> hasBoundaryScvf_;

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

    static constexpr auto dim = Traits::StaticInfo::dim;
    static constexpr auto numScvsPerElement = Traits::StaticInfo::numScvsPerElement;
    static constexpr auto numLateralScvfsPerScv = Traits::StaticInfo::numLateralScvfsPerScv;
    static constexpr auto numLateralScvfsPerElement = Traits::StaticInfo::numLateralScvfsPerElement;
    static constexpr auto minNumScvfsPerElement = Traits::StaticInfo::minNumScvfsPerElement;
    static constexpr auto maxNumScvfsPerElement = Traits::StaticInfo::maxNumScvfsPerElement;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::FCStaggered;
    static constexpr DiscretizationMethod discMethod{};

    static constexpr bool cachingEnabled = false;

    //! export basic grid geometry type for the alternative constructor
    using BasicGridGeometry = BasicGridGeometry_t<GV, Traits>;
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
    //! export the local intersection mapper
    using LocalIntersectionMapper = typename Traits::LocalIntersectionMapper;
    //! export static information
    using StaticInformation = typename Traits::StaticInfo;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;

    //! Constructor with basic grid geometry used to share state with another grid geometry on the same grid view
    FaceCenteredStaggeredFVGridGeometry(std::shared_ptr<BasicGridGeometry> gg, const std::string& paramGroup = "")
    : ParentType(std::move(gg))
    , intersectionMapper_(this->gridView())
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod>::isValid(this->gridView()))
            DUNE_THROW(Dune::InvalidStateException, "The staggered discretization method needs at least an overlap of 1 for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");

        update_();
    }

    //! Constructor from gridView
    FaceCenteredStaggeredFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
    : FaceCenteredStaggeredFVGridGeometry(std::make_shared<BasicGridGeometry>(gridView), paramGroup)
    {}

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

private:

    void update_()
    {
        intersectionMapper_.update(this->gridView());

        // clear local data
        numScvf_ = 0;
        numBoundaryScv_ = 0;
        numBoundaryScvf_ = 0;
        hasBoundaryScvf_.clear();
        scvfIndicesOfElement_.clear();
        outsideVolVarIndices_.clear();

        // determine size of containers
        const auto numElements = this->gridView().size(0);
        scvfIndicesOfElement_.resize(numElements);
        hasBoundaryScvf_.resize(numElements, false);
        numScvs_ = numElements*numScvsPerElement;

        GeometryHelper geometryHelper(this->gridView());

        // get the global scv indices first
        GridIndexType scvfIdx = 0;

        GridIndexType neighborVolVarIdx = numScvs_;

        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);
            assert(numScvsPerElement == element.subEntities(1));

            // the element-wise index sets for finite volume geometry
            auto& globalScvfIndices = scvfIndicesOfElement_[eIdx];
            globalScvfIndices.reserve(maxNumScvfsPerElement);
            globalScvfIndices.resize(minNumScvfsPerElement);

            // keep track of frontal boundary scvfs
            std::size_t numFrontalBoundaryScvfs = 0;

            using LocalIntersectionIndexMapper = FaceCenteredStaggeredLocalIntersectionIndexMapper<GridView>;
            LocalIntersectionIndexMapper localIsMapper;
            localIsMapper.update(this->gridView(), element);

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                const auto localScvIdx = localIsMapper.realToRefIdx(intersection.indexInInside());
                auto localScvfIdx = localScvIdx*(1 + numLateralScvfsPerScv);

                assert(localIsMapper.refToRealIdx(localScvIdx) == intersection.indexInInside());
                // the frontal sub control volume face at the element center
                globalScvfIndices[localScvfIdx] = scvfIdx++;
                ++localScvfIdx;

                if constexpr(dim > 1)
                {
                    // the lateral sub control volume faces
                    for (const auto lateralFacetIndex : Dune::transformedRangeView(geometryHelper.localLaterFaceIndices(localScvIdx),
                                                                                   [&](auto idx) { return localIsMapper.refToRealIdx(idx) ;})
                        )
                    {
                        if (onDomainBoundary_(geometryHelper.intersection(lateralFacetIndex, element)))
                        {
                            outsideVolVarIndices_[scvfIdx] = neighborVolVarIdx++;
                            ++numBoundaryScvf_;
                            hasBoundaryScvf_[eIdx] = true;
                        }

                        globalScvfIndices[localScvfIdx] = scvfIdx++;
                        ++localScvfIdx;
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
            }

            // add global indices of frontal boundary scvfs last
            for (std::size_t i = 0; i < numFrontalBoundaryScvfs; ++i)
                globalScvfIndices.push_back(scvfIdx++);
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

    std::vector<std::vector<GridIndexType>> scvfIndicesOfElement_;

    // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicFaceMap_;
    std::unordered_map<GridIndexType, GridIndexType> outsideVolVarIndices_;
};

} // end namespace Dumux

#endif
