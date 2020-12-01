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
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDamondFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_FV_GRID_GEOMETRY
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_FV_GRID_GEOMETRY

#include <memory>

#include <dumux/common/defaultmappertraits.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/discretization/facecentered/diamond/subcontrolvolume.hh>
#include <dumux/discretization/facecentered/diamond/subcontrolvolumeface.hh>
#include <dumux/discretization/facecentered/diamond/fvelementgeometry.hh>
#include <dumux/discretization/facecentered/diamond/geometryhelper.hh>
#include <dumux/discretization/facecentered/diamond/connectivitymap.hh>


namespace Dumux {

/*!
 * \ingroup XXX
 * \brief The default traits for the xxx finite volume grid geometry
 *        Defines the scv and scvf types and the mapper types
 * \tparam the grid view type TODO docme
 */
template<class GridView>
struct FaceCenteredDiamondDefaultGridGeometryTraits : public DefaultMapperTraits<GridView>
{
    using SubControlVolume = FaceCenteredDiamondSubControlVolume<GridView>;
    using SubControlVolumeFace = FaceCenteredDiamondSubControlVolumeFace<GridView>;
    using IntersectionMapper = ConformingGridIntersectionMapper<GridView>;
    using GeometryHelper = HyperPyramidGeometryHelper<GridView,
                                                      SubControlVolume,
                                                      SubControlVolumeFace>;

    template<class GridGeometry>
    using ConnectivityMap = FaceCenteredDiamondConnectivityMap<GridGeometry>;

    template<class GridGeometry, bool enableCache>
    using LocalView = FaceCenteredDiamondFVElementGeometry<GridGeometry, enableCache>;
};

/*!
 * \ingroup DiamondDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class GridView,
         bool cachingEnabled = false,
         class Traits = FaceCenteredDiamondDefaultGridGeometryTraits<GridView>>
class FaceCenteredDiamondFVGridGeometry;

/*!
 * \ingroup DiamondDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored.
 */
template<class GV, class Traits>
class FaceCenteredDiamondFVGridGeometry<GV, true, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = FaceCenteredDiamondFVGridGeometry<GV, true, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::SmallLocalIndex;
    using Element = typename GV::template Codim<0>::Entity;

    using IntersectionMapper = typename Traits::IntersectionMapper;
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;

    using Scalar = typename GV::ctype;

    static const int dim = GV::dimension;
    static const int dimWorld = GV::dimensionworld;

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::fcstaggered;
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
    FaceCenteredDiamondFVGridGeometry(const GridView& gridView, const std::string& paramGroup = "")
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

    //! update all fvElementGeometries (do this again after grid adaption)
    void update()
    {
        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvIndicesOfElement_.clear();
        scvfIndicesOfElement_.clear();
        intersectionMapper_.update();

        // determine size of containers
        const auto numElements = this->gridView().size(0);
        scvIndicesOfElement_.resize(numElements);
        scvfIndicesOfElement_.resize(numElements);
        numBoundaryScvfOfElement_.resize(numElements, 0);

        outSideBoundaryVolVarIdx_ = 0;
        numBoundaryScv_ = 0;
        numBoundaryScvf_ = 0;

        // determine size of containers
        std::size_t numScvs = 0;
        std::size_t numInteriorScvf = 0;
        GridIndexType scvIdx = 0;
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);
            scvIndicesOfElement_[eIdx].resize(element.subEntities(1));
            numScvs += element.subEntities(1);
            numInteriorScvf += dim==1 ? 2 : 2*element.subEntities(2);

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                scvIndicesOfElement_[eIdx][intersection.indexInInside()] = scvIdx++;
                // handle physical domain boundary
                if (onDomainBoundary_(intersection))
                {
                     ++numBoundaryScv_;
                     ++numBoundaryScvfOfElement_[eIdx];
                }
            }
        }

        numBoundaryScvf_ = numBoundaryScv_;

        std::unordered_map<GridIndexType, GridIndexType> globalScvfIdxToGlobalScvfIdxWithCommonEntity;

         // reserve memory
        scvs_.reserve(numScvs);
        scvfs_.reserve(numInteriorScvf + numBoundaryScvf_);

        // Build the scvs and scv faces
        // get the global scv indices first
        GridIndexType scvfIdx = 0;
        for (const auto& element : elements(this->gridView()))
        {
            GeometryHelper geometryHelper(element);
            const auto eIdx = this->elementMapper().index(element);
            const auto& globalScvIndices = scvIndicesOfElement_[eIdx];

            LocalIndexType localScvfIdx = 0;
            std::vector<GridIndexType> scvfsIndexSet;
            // scvfs are given for each codim 2 entity
            scvfsIndexSet.resize(dim==1 ? 2 : 2*geometryHelper.numInteriorScvf() + numBoundaryScvfOfElement_[eIdx]);

            for (const auto& intersection : intersections(this->gridView(), element))
            {
                const auto dofIndex = intersectionMapper().globalIntersectionIndex(element, intersection.indexInInside());
                const auto intersectionGeometry = intersection.geometry();

                // handle periodic boundaries
                if (onPeriodicBoundary_(intersection))
                {
                    this->setPeriodic();

                    const auto& otherElement = intersection.outside();

                    LocalIndexType otherIntersectionLocalIdx = 0;
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

                if (onDomainBoundary_(intersection))
                {
                    LocalIndexType isIdx = intersection.indexInInside();
                    const auto corners = geometryHelper.getBoundaryScvfCorners(intersection);
                    scvfsIndexSet[localScvfIdx] = scvfIdx++;

                    scvfs_.emplace_back(corners,
                                        intersection,
                                        std::array{isIdx, isIdx},
                                        std::array{globalScvIndices[isIdx], globalScvIndices[isIdx]},
                                        localScvfIdx,
                                        scvfsIndexSet[localScvfIdx]);
                    ++localScvfIdx;
                }

                // the sub control volume
                scvs_.emplace_back(geometryHelper.getScvCorners(intersection.indexInInside()),
                                   intersectionGeometry.center(),
                                   globalScvIndices[intersection.indexInInside()],
                                   intersection.indexInInside(),
                                   dofIndex,
                                   this->elementMapper().index(element),
                                   onDomainBoundary_(intersection));
            }

            for(LocalIndexType idx=0; idx < geometryHelper.numInteriorScvf(); ++idx)
            {
                const auto corners = geometryHelper.getScvfCorners(idx);
                auto scvPair = geometryHelper.getFaceScvPair(idx);
                // the sub control volume faces
                scvfsIndexSet[localScvfIdx] = scvfIdx++;
                auto normal = geometryHelper.normal(scvPair, idx);

                scvfs_.emplace_back(corners,
                                    normal,
                                    std::array{scvPair.first, scvPair.second},
                                    std::array{globalScvIndices[scvPair.first], globalScvIndices[scvPair.second]},
                                    localScvfIdx,
                                    scvfsIndexSet[localScvfIdx]);
                ++localScvfIdx;

                // We duplicate the scvfs because assembly is currently always done with respect to insideScv
                scvfsIndexSet[localScvfIdx] = scvfIdx++;
                scvfs_.emplace_back(corners,
                                    -1*normal,
                                    std::array{scvPair.second, scvPair.first},
                                    std::array{globalScvIndices[scvPair.second], globalScvIndices[scvPair.first]},
                                    localScvfIdx,
                                    scvfsIndexSet[localScvfIdx]);
                ++localScvfIdx;
            }

            // // the frontal sub control volume face at a domain boundary (coincides with element face)
            // if (onDomainBoundary_(intersection))
            // {
            //     ++numBoundaryScvf_;
            //     const auto boundaryCenter = intersectionGeometry.center();
            //     scvfs_.emplace_back(boundaryCenter,
            //                         boundaryCenter,
            //                         std::array{localScvIdx, localOppositeScvIdx}, // TODO outside boundary, periodic, parallel?
            //                         std::array{globalScvIndices[localScvIdx], globalScvIndices[localScvIdx]}, // TODO outside boundary, periodic, parallel?
            //                         localScvfIdx,
            //                         intersectionGeometry.volume(),
            //                         directionIdx,
            //                         sign(intersection.centerUnitOuterNormal()[directionIdx]),
            //                         globalScvfIndices[localScvfIdx],
            //                         0, // should not be used
            //                         SubControlVolumeFace::FaceType::frontal,
            //                         true);
            //     ++localScvfIdx;
            //     hasBoundaryScvf_[eIdx] = true;
            // }
            scvfIndicesOfElement_[eIdx] = std::move(scvfsIndexSet);
        }

        connectivityMap_.update(*this);
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

    //! Get the global sub control volume indices of an element
    const std::vector<GridIndexType>& scvIndicesOfElement(GridIndexType eIdx) const
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
    { return numBoundaryScvfOfElement_[eIdx] > 0; }

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
    std::vector<LocalIndexType> numBoundaryScvfOfElement_;

    std::vector<std::vector<GridIndexType>> scvIndicesOfElement_;
    std::vector<std::vector<GridIndexType>> scvfIndicesOfElement_;

        // a map for periodic boundary vertices
    std::unordered_map<GridIndexType, GridIndexType> periodicFaceMap_;
};

/*!
 * \ingroup DiamondDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored.
 */
template<class GV, class Traits>
class FaceCenteredDiamondFVGridGeometry<GV, false, Traits>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = FaceCenteredDiamondFVGridGeometry<GV, false, Traits>;
    using ParentType = BaseGridGeometry<GV, Traits>;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;
    using Element = typename GV::template Codim<0>::Entity;

    using IntersectionMapper = typename Traits::IntersectionMapper;
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;

public:

};

} // end namespace

#endif
