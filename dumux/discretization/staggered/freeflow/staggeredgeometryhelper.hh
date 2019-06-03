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

namespace Dumux
{

namespace Detail {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Parallel Data stored per sub face
 */
template<class GridView, int upwindSchemeOrder>
struct PairData
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

    std::bitset<upwindSchemeOrder> hasParallelNeighbor;
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
struct AxisData;
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
 * \brief In Axis Data stored per sub face
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
    static constexpr auto dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

    //TODO include assert that checks for quad geometry
    static constexpr auto codimIntersection =  1;
    static constexpr auto codimCommonEntity = 2;
    static constexpr auto numFacetSubEntities = (dim == 2) ? 2 : 4;
    static constexpr auto numfacets = dimWorld * 2;
    static constexpr auto numPairs = 2 * (dimWorld - 1);

    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

public:
    using PairData = Detail::PairData<GridView, upwindSchemeOrder>;
    using AxisData = Detail::AxisData<GridView, upwindSchemeOrder>;

    FreeFlowStaggeredGeometryHelper(const Element& element, const GridView& gridView) : element_(element), elementGeometry_(element.geometry()), gridView_(gridView)
    { }

    //! update the local face
    template<class IntersectionMapper>
    void updateLocalFace(const IntersectionMapper& intersectionMapper_, const Intersection& intersection)
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
     * \brief Returns the dirction index of the primary facet (0 = x, 1 = y, 2 = z)
     */
    unsigned int directionIndex() const
    {
        return Dumux::directionIndex(intersection_.centerUnitOuterNormal());
    }

    /*!
     * \brief Returns the dirction index of the facet passed as an argument (0 = x, 1 = y, 2 = z)
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
        fillOuterAxisData_(std::integral_constant<bool, useHigherOrder>{});

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
     * \brief Fills the axis data with the outer Dofs and Distances.
     *        For first order upwind methods no outer information is required.
     */
    void fillOuterAxisData_(std::false_type) {}

    /*!
     * \brief Fills the axis data with the outer dofs and distances.
     *        This function builds and extended stencil and is therefore only called when higher order upwinding methods are prescribed.
     */
    void fillOuterAxisData_(std::true_type)
    {
        // reset the axis data struct
        axisData_ = {};

        const auto numForwardBackwardAxisDofs = axisData_.inAxisForwardDofs.size();

        // Set the Forward Dofs
        std::stack<Element> inAxisForwardElementStack;
        const auto inIdx = intersection_.indexInInside();
        auto selfFace = getFacet_(inIdx, element_);
        if(intersection_.neighbor())
        {
            inAxisForwardElementStack.push(intersection_.outside());
            bool keepStackingForward = (inAxisForwardElementStack.size() < numForwardBackwardAxisDofs);
            while(keepStackingForward)
            {
                auto e = inAxisForwardElementStack.top();
                for(const auto& intersection : intersections(gridView_,e))
                {
                    if( (intersection.indexInInside() == inIdx ) )
                    {
                        if( intersection.neighbor())
                        {
                            inAxisForwardElementStack.push(intersection.outside());
                            keepStackingForward = (inAxisForwardElementStack.size() < numForwardBackwardAxisDofs);
                        }
                        else
                        {
                            keepStackingForward = false;
                        }
                    }
                }
            }
        }

        std::vector<GlobalPosition> forwardFaceCoordinates(inAxisForwardElementStack.size(), selfFace.geometry().center());
        while(!inAxisForwardElementStack.empty())
        {
            int forwardIdx = inAxisForwardElementStack.size()-1;
            axisData_.hasForwardNeighbor.set(forwardIdx, true);
            axisData_.inAxisForwardDofs[forwardIdx] = gridView_.indexSet().subIndex(inAxisForwardElementStack.top(), inIdx, codimIntersection);
            selfFace = getFacet_(inIdx, inAxisForwardElementStack.top());
            forwardFaceCoordinates[forwardIdx] = selfFace.geometry().center();
            inAxisForwardElementStack.pop();
        }

        const auto self = getFacet_(inIdx, element_);
        for (int i = 0; i< forwardFaceCoordinates.size(); i++)
        {
            if (i == 0)
                axisData_.inAxisForwardDistances[i] = (forwardFaceCoordinates[i] - self.geometry().center()).two_norm();
            else
                axisData_.inAxisForwardDistances[i] = (forwardFaceCoordinates[i]- forwardFaceCoordinates[i-1]).two_norm();
        }

        // Set the Backward Dofs
        std::stack<Element> inAxisBackwardElementStack;
        const auto oppIdx = localOppositeIdx_(inIdx);
        auto oppFace = getFacet_(oppIdx, element_);
        for(const auto& intersection : intersections(gridView_, element_))
        {
            if(intersection.indexInInside() == oppIdx)
            {
                if(intersection.neighbor())
                {
                    inAxisBackwardElementStack.push(intersection.outside());
                    bool keepStackingBackward = (inAxisBackwardElementStack.size() < numForwardBackwardAxisDofs);
                    while(keepStackingBackward)
                    {
                        auto e = inAxisBackwardElementStack.top();
                        for(const auto& intersectionOut : intersections(gridView_,e))
                        {

                            if( (intersectionOut.indexInInside() == oppIdx ) )
                            {
                                if( intersectionOut.neighbor())
                                {
                                    inAxisBackwardElementStack.push(intersectionOut.outside());
                                    keepStackingBackward = (inAxisBackwardElementStack.size() < numForwardBackwardAxisDofs);
                                }
                                else
                                {
                                    keepStackingBackward = false;
                                }
                            }
                        }
                    }
                }
            }
        }

        std::vector<GlobalPosition> backwardFaceCoordinates(inAxisBackwardElementStack.size(), oppFace.geometry().center());
        while(!inAxisBackwardElementStack.empty())
        {
            int backwardIdx = inAxisBackwardElementStack.size()-1;
            axisData_.hasBackwardNeighbor.set(backwardIdx, true);
            axisData_.inAxisBackwardDofs[backwardIdx] = gridView_.indexSet().subIndex(inAxisBackwardElementStack.top(), oppIdx, codimIntersection);
            oppFace = getFacet_(oppIdx, inAxisBackwardElementStack.top());
            backwardFaceCoordinates[backwardIdx] = oppFace.geometry().center();
            inAxisBackwardElementStack.pop();
        }

        const auto opposite = getFacet_(oppIdx, element_);
        for (int i = 0; i< backwardFaceCoordinates.size(); i++)
        {
            if (i == 0)
                axisData_.inAxisBackwardDistances[i] = (backwardFaceCoordinates[i] - opposite.geometry().center()).two_norm();
            else
                axisData_.inAxisBackwardDistances[i] = (backwardFaceCoordinates[i] - backwardFaceCoordinates[i-1]).two_norm();
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
                pairData_[numPairInnerLateralIdx].lateralStaggeredFaceCenter = std::move(selfFacetCenter + distance);

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
                    pairData_[numPairOuterLateralIdx].lateralDistance = std::move(normalDistanceoffset.two_norm());
                    numPairOuterLateralIdx++;
                }
            }
        }

        fillParallelPairData_(std::integral_constant<bool, useHigherOrder>{});
    }

    /*!
     * \brief Fills the pair data with the parallel dofs and distances
     *        This function is only called when the simple first order upwind methods are used.
     */
    void fillParallelPairData_(std::false_type)
    {
        // set basic global positions and stencil size definitions
        // get the parallel Dofs
        const auto parallelLocalIdx = intersection_.indexInInside();
        SmallLocalIndexType numPairParallelIdx = 0;

        for (const auto& intersection : intersections(gridView_, element_))
        {
            if (facetIsNormal_(intersection.indexInInside(), parallelLocalIdx))
            {
                auto parallelAxisIdx = directionIndex(intersection);
                if (intersection.neighbor())
                {
                    // If the lateral intersection has a neighboring cell, go in and store the parallel information.
                    const auto& outerElement = intersection.outside();
                    pairData_[numPairParallelIdx].hasParallelNeighbor.set(0, true);
                    pairData_[numPairParallelIdx].parallelDofs[0] = gridView_.indexSet().subIndex(outerElement, parallelLocalIdx, codimIntersection);
                    pairData_[numPairParallelIdx].parallelCellWidths[0] = setParallelPairCellWidths_(outerElement, parallelAxisIdx);
                }

                numPairParallelIdx++;
            }
        }
    }

    /*!
     * \brief Fills the pair data with the parallel dofs and distances.
     *        This function builds and extended stencil and is therefore only called when higher order upwinding methods are prescribed.
     */
    void fillParallelPairData_(std::true_type)
    {
        // set basic global positions and stencil size definitions
        const auto numParallelFaces = pairData_[0].parallelCellWidths.size();

        // get the parallel Dofs
        const auto parallelLocalIdx = intersection_.indexInInside();
        SmallLocalIndexType numPairParallelIdx = 0;
        std::stack<Element> parallelElementStack;
        for(const auto& intersection : intersections(gridView_, element_))
        {
            if( facetIsNormal_(intersection.indexInInside(), parallelLocalIdx) )
            {
                if( intersection.neighbor() )
                {
                    auto parallelAxisIdx = directionIndex(intersection);
                    auto localLateralIntersectionIndex = intersection.indexInInside();
                    auto e = element_;

                    bool keepStacking =  (parallelElementStack.size() < numParallelFaces);
                    while(keepStacking)
                    {
                        for(const auto& lateralIntersection : intersections(gridView_, e))
                        {
                            if( lateralIntersection.indexInInside() == localLateralIntersectionIndex )
                            {
                                if( lateralIntersection.neighbor() )
                                {
                                    parallelElementStack.push(lateralIntersection.outside());
                                    keepStacking = (parallelElementStack.size() < numParallelFaces);
                                }
                                else
                                {
                                    keepStacking = false;
                                }
                            }
                        }
                        e = parallelElementStack.top();
                    }

                    while(!parallelElementStack.empty())
                    {
                        pairData_[numPairParallelIdx].hasParallelNeighbor.set(parallelElementStack.size()-1, true);
                        pairData_[numPairParallelIdx].parallelDofs[parallelElementStack.size()-1] = gridView_.indexSet().subIndex(parallelElementStack.top(), parallelLocalIdx, codimIntersection);
                        pairData_[numPairParallelIdx].parallelCellWidths[parallelElementStack.size()-1] = setParallelPairCellWidths_(parallelElementStack.top(), parallelAxisIdx);
                        parallelElementStack.pop();
                    }

                }

                numPairParallelIdx++;
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
    };

    auto getFacet_(const int localFacetIdx, const Element& element) const
    {
        return element.template subEntity <1> (localFacetIdx);
    };

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

        // store the element distance
        const auto& outerLateralFacet = getFacet_(isIdx, element);
        const auto outerLateralFacetPos = outerLateralFacet.geometry().center();
        const auto& innerLateralFacet = getFacet_(isIdx, element_);
        const auto innerLateralFacetPos = innerLateralFacet.geometry().center();
        pairData_[numPairsIdx].lateralDistance = (innerLateralFacetPos - outerLateralFacetPos).two_norm();
    }

    //! Sets the information about the parallel distances
    Scalar setParallelPairCellWidths_(const Element& element, const int parallelAxisIdx) const
    {
        std::vector<GlobalPosition> faces;
        faces.reserve(numfacets);
        for (const auto& intElement : intersections(gridView_, element))
        {
            faces.push_back(intElement.geometry().center());
        }
        switch(parallelAxisIdx)
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
    const Element element_; //!< The respective element
    const typename Element::Geometry elementGeometry_; //!< Reference to the element geometry
    const GridView gridView_; //!< The grid view
    AxisData axisData_; //!< Data related to forward and backward faces
    std::array<PairData, numPairs> pairData_; //!< Collection of pair information related to lateral and parallel faces
};

} // end namespace Dumux

#endif
