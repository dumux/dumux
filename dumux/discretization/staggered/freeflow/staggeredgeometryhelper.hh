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
 * \copydoc Dumux::FreeFlowStaggeredGeometryHelper
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GEOMETRY_HELPER_HH

#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>
#include <type_traits>
#include <algorithm>
#include <array>
#include <bitset>

namespace Dumux {
namespace Detail {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Parallel Data stored per sub face
 *
 * ------------
 * |          |
 * |          |
 * |          |
 * -----------------------
 * | yyyyyyyy s          |
 * | yyyyyyyy s          |
 * | yyyyyyyy s          |
 * -----------------------
 * In this corner geometry, hasParallelNeighbor will return true for subcontrolvolumeface s belonging to the
 * element filled by 'y's, but hasParallelNeighbor will return false for the subcontrolvolumeface that has the
 * same dofIndex. We name this situation hasHalfParallelNeighbor.
 *
 * ------------
 * | yyyyyyyy s
 * | yyyyyyyy s
 * | yyyyyyyy s
 * -----------------------
 * |          |          |
 * |          |          |
 * |          |          |
 * -----------------------
 * In this corner geometry, hasParallelNeighbor will return true for subcontrolvolumeface s belonging to the
 * element filled by 'y's. However, as there also might be a boundary velocity value known at the corner, which
 * can be used instead of the standard parallel velocity in some cases, we want to identify this situation. We
 * name it cornerParallelNeighbor.
 */
template<class GridView, int upwindSchemeOrder>
struct PairData
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

    std::bitset<upwindSchemeOrder> hasParallelNeighbor;
    bool hasHalfParallelNeighbor = false;
    bool hasCornerParallelNeighbor = false;
    std::array<GridIndexType, upwindSchemeOrder> parallelDofs;
    std::array<Scalar, upwindSchemeOrder> parallelCellWidths;
    bool hasOuterLateral = false;
    std::pair<GridIndexType, GridIndexType> lateralPair;
    SmallLocalIndexType localLateralFaceIdx;
    Scalar lateralDistance;
    GlobalPosition lateralStaggeredFaceCenter;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief In Axis Data stored per sub face
 */
template<class GridView, int upwindSchemeOrder>
struct AxisData
{
    using Scalar = typename GridView::ctype;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    GridIndexType selfDof;
    GridIndexType oppositeDof;
    std::bitset<upwindSchemeOrder-1> hasForwardNeighbor;
    std::bitset<upwindSchemeOrder-1> hasBackwardNeighbor;
    std::array<GridIndexType, upwindSchemeOrder-1> inAxisForwardDofs;
    std::array<GridIndexType, upwindSchemeOrder-1> inAxisBackwardDofs;
    Scalar selfToOppositeDistance;
    std::array<Scalar, upwindSchemeOrder-1> inAxisForwardDistances;
    std::array<Scalar, upwindSchemeOrder-1> inAxisBackwardDistances;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief In Axis Data stored per sub face for first-order scheme
 */
template<class GridView>
struct AxisData<GridView, 1>
{
    using Scalar = typename GridView::ctype;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    GridIndexType selfDof;
    GridIndexType oppositeDof;
    Scalar selfToOppositeDistance;
};

} // namespace Detail

/*!
 * \ingroup StaggeredDiscretization
 * \brief Returns the dirction index of the facet (0 = x, 1 = y, 2 = z)
 */
template<class Vector>
inline static unsigned int directionIndex(Vector&& vector)
{
    const auto eps = 1e-8;
    return std::find_if(vector.begin(), vector.end(), [eps](const auto& x) { return std::abs(x) > eps; } ) - vector.begin();
}

/*!
 * \ingroup StaggeredDiscretization
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the free flow staggered discretization method.
 */
template<class GridView, int upwindSchemeOrder>
class FreeFlowStaggeredGeometryHelper
{
    using Scalar = typename GridView::ctype;
    static constexpr auto dim = GridView::dimension;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

    //TODO include assert that checks for quad geometry
    static constexpr auto codimIntersection =  1;
    static constexpr auto numFacets = dim * 2;
    static constexpr auto numPairs = 2 * (dim - 1);

    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

public:
    using PairData = Detail::PairData<GridView, upwindSchemeOrder>;
    using AxisData = Detail::AxisData<GridView, upwindSchemeOrder>;

    FreeFlowStaggeredGeometryHelper(const Element& element, const GridView& gridView)
    : element_(element)
    , gridView_(gridView)
    { }

    //! update the local face
    template<class IntersectionMapper>
    void updateLocalFace(const IntersectionMapper&, const Intersection& intersection)
    {
        intersection_ = intersection;
        fillAxisData_();
        fillPairData_();
    }

    /*!
     * \brief Returns the local index of the face (i.e. the intersection)
     */
    SmallLocalIndexType localFaceIndex() const
    {
        return intersection_.indexInInside();
    }

    /*!
     * \brief Returns a copy of the axis data
     */
    AxisData axisData() const
    {
        return axisData_;
    }

    /*!
     * \brief Returns a copy of the pair data
     */
    std::array<PairData, numPairs> pairData() const
    {
        return pairData_;
    }

    /*!
     * \brief Returns the direction index of the primary facet (0 = x, 1 = y, 2 = z)
     */
    unsigned int directionIndex() const
    {
        return Dumux::directionIndex(intersection_.centerUnitOuterNormal());
    }

    /*!
     * \brief Returns the direction index of the facet passed as an argument (0 = x, 1 = y, 2 = z)
     */
    unsigned int directionIndex(const Intersection& intersection) const
    {
        return Dumux::directionIndex(std::move(intersection.centerUnitOuterNormal()));
    }

private:

    /*!
     * \brief Fills all entries of the in axis data
     *        Calls a function to extend the axis data for higher order upwind methods if required.
     */
    void fillAxisData_()
    {
        if constexpr (useHigherOrder)
            fillOuterAxisData_();

        const auto inIdx = intersection_.indexInInside();
        const auto oppIdx = localOppositeIdx_(inIdx);

        // Set the self Dof
        axisData_.selfDof = gridView_.indexSet().subIndex(intersection_.inside(), inIdx, codimIntersection);

        // Set the opposite Dof
        axisData_.oppositeDof = gridView_.indexSet().subIndex(intersection_.inside(), oppIdx, codimIntersection);

        // Set the Self to Opposite Distance
        const auto self = getFacet_(inIdx, element_);
        const auto opposite = getFacet_(oppIdx, element_);
        axisData_.selfToOppositeDistance = (self.geometry().center() - opposite.geometry().center()).two_norm();

    }

    /*!
     * \brief Fills the axis data with the outer dofs and distances.
     *        This function builds and extended stencil (see AxisData for upwindSchemeOrder > 1)
     *        and is therefore only called when higher order upwinding methods are prescribed.
     */
    void fillOuterAxisData_()
    {
        // reset the axis data struct
        axisData_ = {};

        // Set the forward dofs
        const auto thisFaceLocalIdx = intersection_.indexInInside();
        const auto& faceCenterPos = getFacet_(thisFaceLocalIdx, element_).geometry().center();
        if (intersection_.neighbor())
            addForwardNeighborAxisData_(intersection_.outside(), 0, thisFaceLocalIdx, faceCenterPos);

        // Set the backward dofs
        const auto oppositeFaceLocalIdx = localOppositeIdx_(thisFaceLocalIdx);
        const auto& oppositeFaceCenterPos = getFacet_(oppositeFaceLocalIdx, element_).geometry().center();
        for (const auto& intersection : intersections(gridView_, element_))
        {
            if (intersection.indexInInside() == oppositeFaceLocalIdx && intersection.neighbor())
            {
                addBackwardNeighborAxisData_(intersection.outside(), 0, oppositeFaceLocalIdx, oppositeFaceCenterPos);
                break;
            }
        }
    }

    /*!
     * \brief Recursively add axis data in forward direction
     */
    void addForwardNeighborAxisData_(const Element& neighbor,
                                     const std::size_t distance,
                                     const unsigned int localIntersectionIndex,
                                     const GlobalPosition& faceCenterPos)
    {
        // fill the forward axis data for this element in the stencil
        axisData_.hasForwardNeighbor.set(distance, true);
        axisData_.inAxisForwardDofs[distance] = gridView_.indexSet().subIndex(neighbor, localIntersectionIndex, codimIntersection);
        const auto& forwardFacePos = getFacet_(localIntersectionIndex, neighbor).geometry().center();
        axisData_.inAxisForwardDistances[distance] = (forwardFacePos - faceCenterPos).two_norm();

        // recursion end if we added all axis data for the given stencil size (numParallelFaces)
        const auto numForwardBackwardAxisDofs = axisData_.inAxisForwardDofs.size();
        if (distance >= numForwardBackwardAxisDofs-1)
            return;

        for (const auto& intersection : intersections(gridView_, neighbor))
        {
            if (intersection.indexInInside() == localIntersectionIndex && intersection.neighbor())
            {
                addForwardNeighborAxisData_(intersection.outside(), distance+1, localIntersectionIndex, forwardFacePos);
                break;
            }
        }
    }

    /*!
     * \brief Recursively add axis data in backward direction
     */
    void addBackwardNeighborAxisData_(const Element& neighbor,
                                      const std::size_t distance,
                                      const unsigned int localIntersectionIndex,
                                      const GlobalPosition& faceCenterPos)
    {
        // fill the forward axis data for this element in the stencil
        axisData_.hasBackwardNeighbor.set(distance, true);
        axisData_.inAxisBackwardDofs[distance] = gridView_.indexSet().subIndex(neighbor, localIntersectionIndex, codimIntersection);
        const auto& backwardFacePos = getFacet_(localIntersectionIndex, neighbor).geometry().center();
        axisData_.inAxisBackwardDistances[distance] = (backwardFacePos - faceCenterPos).two_norm();

        // recursion end if we added all axis data for the given stencil size (numParallelFaces)
        const auto numForwardBackwardAxisDofs = axisData_.inAxisForwardDofs.size();
        if (distance >= numForwardBackwardAxisDofs-1)
            return;

        for (const auto& intersection : intersections(gridView_, neighbor))
        {
            if (intersection.indexInInside() == localIntersectionIndex && intersection.neighbor())
            {
                addBackwardNeighborAxisData_(intersection.outside(), distance+1, localIntersectionIndex, backwardFacePos);
                break;
            }
        }
    }

    /*!
     * \brief Fills the pair data with the lateral dofs and distances
     *        and calls a further function to collect the parallel dofs and distances
     */
    void fillPairData_()
    {
        // reset the pair data structs
        pairData_ = {};

        // set basic global positions
        const auto& selfElementCenter = element_.geometry().center();
        const auto& selfFacetCenter = intersection_.geometry().center();

        // get the inner lateral Dof Index
        SmallLocalIndexType numPairInnerLateralIdx = 0;
        for (const auto& innerElementIntersection : intersections(gridView_, element_))
        {
            if (facetIsNormal_(innerElementIntersection.indexInInside(), intersection_.indexInInside()))
            {
                const auto innerElementIntersectionIdx = innerElementIntersection.indexInInside();
                setLateralPairFirstInfo_(innerElementIntersectionIdx, element_, numPairInnerLateralIdx);

                const auto distance = innerElementIntersection.geometry().center() - selfElementCenter;
                pairData_[numPairInnerLateralIdx].lateralStaggeredFaceCenter = selfFacetCenter+distance;

                numPairInnerLateralIdx++;
            }
        }

        // get the outer lateral Dof Index
        SmallLocalIndexType numPairOuterLateralIdx = 0;
        if (intersection_.neighbor())
        {
            // the direct neighbor element and the respective intersection index
            const auto& directNeighborElement = intersection_.outside();

            for (const auto& directNeighborElementIntersection : intersections(gridView_, directNeighborElement))
            {
                // skip the directly neighboring face itself and its opposing one
                if (facetIsNormal_(directNeighborElementIntersection.indexInInside(), intersection_.indexInOutside()))
                {
                    const auto directNeighborElemIsIdx = directNeighborElementIntersection.indexInInside();
                    setLateralPairSecondInfo_(directNeighborElemIsIdx, directNeighborElement, numPairOuterLateralIdx);
                    numPairOuterLateralIdx++;
                }
            }
        }
        else // intersection is on boundary
        {
            for (const auto& intersection : intersections(gridView_, element_))
            {
                if (facetIsNormal_(intersection.indexInInside(), intersection_.indexInInside()))
                {
                    assert(!pairData_[numPairOuterLateralIdx].hasOuterLateral);

                    const auto normalDistanceoffset = selfFacetCenter - selfElementCenter;
                    pairData_[numPairOuterLateralIdx].lateralDistance = normalDistanceoffset.two_norm();
                    numPairOuterLateralIdx++;
                }
            }
        }

        // fill the pair data with the parallel dofs and distances
        const auto parallelLocalIdx = intersection_.indexInInside();
        SmallLocalIndexType numPairParallelIdx = 0;
        for (const auto& intersection : intersections(gridView_, element_))
        {
            if (facetIsNormal_(intersection.indexInInside(), parallelLocalIdx))
            {
                if (intersection.neighbor())
                {
                    // recursively insert parallel neighbor faces into pair data
                    const auto parallelAxisIdx = directionIndex(intersection);
                    const auto localLateralIntersectionIndex = intersection.indexInInside();

                   /*
                    *       ------------
                    *       |          |
                    *       |          |
                    *       |          |
                    *       iiiiiiiiiii*bbbbbbbbbbb
                    *       |          o zzzzzzzz |
                    *       |          o zzzzzzzz |
                    *       |          o zzzzzzzz |
                    *       -----------------------
                    *
                    *       i:intersection,o:intersection_, b: outerIntersection, z: intersection_.outside()
                    */
                   if (intersection_.neighbor())
                        for (const auto& outerIntersection : intersections(gridView_, intersection_.outside()))
                            if (intersection.indexInInside() == outerIntersection.indexInInside())
                                if (!outerIntersection.neighbor())
                                    pairData_[numPairParallelIdx].hasHalfParallelNeighbor = true;

                   /*       ------------
                    *       |          o
                    *       |          o
                    *       |          o
                    *       iiiiiiiiiii------------
                    *       | zzzzzzzz b          |
                    *       | zzzzzzzz b          |
                    *       | zzzzzzzz b          |
                    *       -----------------------
                    *
                    *       i:intersection,o:intersection_, b: outerIntersection, z: intersection.outside()
                    */
                   if (!intersection_.neighbor())
                        for (const auto& outerIntersection : intersections(gridView_, intersection.outside()))
                            if (intersection_.indexInInside() == outerIntersection.indexInInside())
                                if (outerIntersection.neighbor())
                                    pairData_[numPairParallelIdx].hasCornerParallelNeighbor = true;


                    addParallelNeighborPairData_(intersection.outside(), 0, localLateralIntersectionIndex, parallelLocalIdx, parallelAxisIdx, numPairParallelIdx);
                }

                ++numPairParallelIdx;
            }
        }
    }

    // zero distance means direct neighbor to element_
    void addParallelNeighborPairData_(const Element& neighbor,
                                      const std::size_t distance,
                                      const unsigned int localLateralIntersectionIndex,
                                      const unsigned int parallelLocalIdx,
                                      const unsigned int parallelAxisIdx,
                                      const SmallLocalIndexType dataIdx)
    {
        // fill the pair data for this element in the stencil
        pairData_[dataIdx].hasParallelNeighbor.set(distance, true);
        pairData_[dataIdx].parallelDofs[distance] = gridView_.indexSet().subIndex(neighbor, parallelLocalIdx, codimIntersection);
        pairData_[dataIdx].parallelCellWidths[distance] = setParallelPairCellWidths_(neighbor, parallelAxisIdx);

        // recursion end if we added all parallel data for the given stencil size (numParallelFaces)
        const auto numParallelFaces = pairData_[0].parallelCellWidths.size();
        if (distance >= numParallelFaces-1)
            return;

        // continue with neighbor's neighbor element in the same direction, if we find one
        for (const auto& lateralIntersection : intersections(gridView_, neighbor))
        {
            // we only expect to find one intersection matching this condition
            if (lateralIntersection.indexInInside() == localLateralIntersectionIndex)
            {
                // if there is no neighbor recursion ends here
                if (lateralIntersection.neighbor())
                    addParallelNeighborPairData_(lateralIntersection.outside(), distance+1,
                                                 localLateralIntersectionIndex, parallelLocalIdx, parallelAxisIdx, dataIdx);
                break;
            }
        }
    }

    /*!
     * \brief Returns the local opposing intersection index
     *
     * \param idx The local index of the intersection itself
     */
    int localOppositeIdx_(const int idx) const
    {
        return (idx % 2) ? (idx - 1) : (idx + 1);
    }

    /*!
     * \brief Returns true if the intersection lies normal to another given intersection
     *
     * \param selfIdx The local index of the intersection itself
     * \param otherIdx The local index of the other intersection
     */
    bool facetIsNormal_(const int selfIdx, const int otherIdx) const
    {
        return !(selfIdx == otherIdx || localOppositeIdx_(selfIdx) == otherIdx);
    }

    auto getFacet_(const int localFacetIdx, const Element& element) const
    {
        return element.template subEntity <1> (localFacetIdx);
    }

    //! Sets the information about the lateral faces (within the element)
    void setLateralPairFirstInfo_(const int isIdx, const Element& element, const int numPairsIdx)
    {
        // store the inner lateral dofIdx
        pairData_[numPairsIdx].lateralPair.first = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);

        // store the local lateral facet index
        pairData_[numPairsIdx].localLateralFaceIdx = isIdx;
    }

    //! Sets the information about the lateral faces (in the neighbor element)
    void setLateralPairSecondInfo_(const int isIdx, const Element& element, const int numPairsIdx)
    {
        // store the dofIdx
        pairData_[numPairsIdx].hasOuterLateral = true;
        pairData_[numPairsIdx].lateralPair.second = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);

        // set basic global positions
        const auto& selfFacetCenter = intersection_.geometry().center();
        const auto& selfElementCenter = element_.geometry().center();

        const auto& neighborElement = intersection_.outside();
        const auto& neighborElementCenter = neighborElement.geometry().center();
        const auto& neighborFacetCenter = getFacet_(intersection_.indexInOutside(), neighborElement).geometry().center();

        const Scalar insideLateralDistance = (selfFacetCenter - selfElementCenter).two_norm();
        const Scalar outsideLateralDistance = (neighborFacetCenter - neighborElementCenter).two_norm();

        pairData_[numPairsIdx].lateralDistance = insideLateralDistance + outsideLateralDistance;
    }

    //! Sets the information about the parallel distances
    Scalar setParallelPairCellWidths_(const Element& element, const int parallelAxisIdx) const
    {
        std::vector<GlobalPosition> faces;
        faces.reserve(numFacets);
        for (const auto& intersection : intersections(gridView_, element))
            faces.push_back(intersection.geometry().center());

        switch (parallelAxisIdx)
        {
            case 0:
                return (faces[1] - faces[0]).two_norm();
            case 1:
                return (faces[3] - faces[2]).two_norm();
            case 2:
            {
                assert(dim == 3);
                return (faces[5] - faces[4]).two_norm();
            }
            default:
                DUNE_THROW(Dune::InvalidStateException, "Something went terribly wrong");
        }
    }

    Intersection intersection_; //!< The intersection of interest
    const Element& element_; //!< The respective element
    const GridView gridView_; //!< The grid view
    AxisData axisData_; //!< Data related to forward and backward faces
    std::array<PairData, numPairs> pairData_; //!< Collection of pair information related to lateral and parallel faces
};

} // end namespace Dumux

#endif
