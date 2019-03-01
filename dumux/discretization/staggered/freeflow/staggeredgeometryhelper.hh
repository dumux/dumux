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
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/math.hh>
#include <type_traits>
#include <algorithm>
#include <array>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief Parallel Data stored per sub face
 */
template<class Scalar, class GlobalPosition, int upwindSchemeOrder>
struct PairData
{
    std::array<int, upwindSchemeOrder> parallelDofs;
    std::array<Scalar, upwindSchemeOrder+1> parallelDistances; // TODO store only two distances, not three.
    std::pair<signed int, signed int> normalPair;
    int localNormalFaceIdx;
    Scalar normalDistance;
    GlobalPosition virtualFirstParallelFaceDofPos;
};


/*!
 * \ingroup StaggeredDiscretization
 * \brief In Axis Data stored per sub face
 */
template<class Scalar, int upwindSchemeOrder>
struct AxisData;
/*!
 * \ingroup StaggeredDiscretization
 * \brief In Axis Data stored per sub face
 */
template<class Scalar, int upwindSchemeOrder>
struct AxisData
{
    AxisData()
    {
        inAxisForwardDofs.fill(-1);
        inAxisBackwardDofs.fill(-1);
        inAxisForwardDistances.fill(0);
        inAxisBackwardDistances.fill(0);
    }

    int selfDof;
    int oppositeDof;
    std::array<int, upwindSchemeOrder-1> inAxisForwardDofs;
    std::array<int, upwindSchemeOrder-1> inAxisBackwardDofs;
    Scalar selfToOppositeDistance;
    std::array<Scalar, upwindSchemeOrder-1> inAxisForwardDistances;
    std::array<Scalar, upwindSchemeOrder-1> inAxisBackwardDistances;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief In Axis Data stored per sub face
 */
template<class Scalar>
struct AxisData<Scalar, 1>
{
    int selfDof;
    int oppositeDof;
    Scalar selfToOppositeDistance;
};

// template<class Scalar, int upwindSchemeOrder>
// using AxisData

/*!
 * \ingroup StaggeredDiscretization
 * \brief Returns the dirction index of the facet (0 = x, 1 = y, 2 = z)
 */
template<class Vector>
inline static unsigned int directionIndex(Vector&& vector)
{
    const auto eps = 1e-8;
    const int idx = std::find_if(vector.begin(), vector.end(), [eps](const auto& x) { return std::abs(x) > eps; } ) - vector.begin();
    return idx;
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
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

    //TODO include assert that checks for quad geometry
    static constexpr int codimIntersection =  1;
    static constexpr int codimCommonEntity = 2;
    static constexpr int numFacetSubEntities = (dim == 2) ? 2 : 4;
    static constexpr int numfacets = dimWorld * 2;
    static constexpr int numPairs = 2 * (dimWorld - 1);

    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;


public:

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
    int localFaceIndex() const
    {
        return intersection_.indexInInside();
    }

    /*!
     * \brief Returns a copy of the axis data
     */
    auto axisData() const
    {
        return axisData_;
    }

    /*!
     * \brief Returns a copy of the pair data
     */
    auto pairData() const
    {
        return pairData_;
    }

    /*!
     * \brief Returns the dirction index of the primary facet (0 = x, 1 = y, 2 = z)
     */
    int directionIndex() const
    {
        return Dumux::directionIndex(std::move(intersection_.centerUnitOuterNormal()));
    }

    /*!
     * \brief Returns the dirction index of the facet passed as an argument (0 = x, 1 = y, 2 = z)
     */
    int directionIndex(const Intersection& intersection) const
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

        fillOuterAxisData_(std::integral_constant<bool, useHigherOrder>{});
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
            axisData_.inAxisForwardDofs[forwardIdx] = gridView_.indexSet().subIndex(inAxisForwardElementStack.top(), inIdx, codimIntersection);
            selfFace = getFacet_(inIdx, inAxisForwardElementStack.top());
            forwardFaceCoordinates[forwardIdx] = selfFace.geometry().center();
            inAxisForwardElementStack.pop();
        }

        const auto self = getFacet_(inIdx, element_);
        for(int i = 0; i< forwardFaceCoordinates.size(); i++)
        {
            if(i == 0)
            {
                axisData_.inAxisForwardDistances[i] = (forwardFaceCoordinates[i] - self.geometry().center()).two_norm();
            }
            else
            {
                axisData_.inAxisForwardDistances[i] = (forwardFaceCoordinates[i]- forwardFaceCoordinates[i-1]).two_norm();
            }
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
            this->axisData_.inAxisBackwardDofs[backwardIdx] = gridView_.indexSet().subIndex(inAxisBackwardElementStack.top(), oppIdx, codimIntersection);
            oppFace = getFacet_(oppIdx, inAxisBackwardElementStack.top());
            backwardFaceCoordinates[backwardIdx] = oppFace.geometry().center();
            inAxisBackwardElementStack.pop();
        }

        const auto opposite = getFacet_(oppIdx, element_);
        for(int i = 0; i< backwardFaceCoordinates.size(); i++)
        {
            if(i == 0)
            {
                this->axisData_.inAxisBackwardDistances[i] = (backwardFaceCoordinates[i] - opposite.geometry().center()).two_norm();
            }
            else
            {
                this->axisData_.inAxisBackwardDistances[i] = (backwardFaceCoordinates[i] - backwardFaceCoordinates[i-1]).two_norm();
            }
        }
    }

    /*!
     * \brief Fills the pair data with the normal dofs and distances
     *        and calls a further function to collect the parallel dofs and distances
     */
    void fillPairData_()
    {
        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(auto& data : pairData_)
        {
            // outer normals
            data.normalPair.second = -1;
        }

        // set basic global positions
        const auto& elementCenter = this->element_.geometry().center();
        const auto& FacetCenter = intersection_.geometry().center();

        // get the inner normal Dof Index
        int numPairInnerNormalIdx = 0;
        for(const auto& innerElementIntersection : intersections(gridView_, element_))
        {
            if(facetIsNormal_(innerElementIntersection.indexInInside(), intersection_.indexInInside()))
            {
                const int innerElementIntersectionIdx = innerElementIntersection.indexInInside();
                setNormalPairFirstInfo_(innerElementIntersectionIdx, element_, numPairInnerNormalIdx);
                numPairInnerNormalIdx++;
            }
        }

        // get the outer normal Dof Index
        int numPairOuterNormalIdx = 0;
        if(intersection_.neighbor())
        {
            // the direct neighbor element and the respective intersection index
            const auto& directNeighborElement = intersection_.outside();

            for(const auto& directNeighborElementIntersection : intersections(gridView_, directNeighborElement))
            {
                // skip the directly neighboring face itself and its opposing one
                if(facetIsNormal_(directNeighborElementIntersection.indexInInside(), intersection_.indexInOutside()))
                {
                    const int directNeighborElemIsIdx = directNeighborElementIntersection.indexInInside();
                    setNormalPairSecondInfo_(directNeighborElemIsIdx, directNeighborElement, numPairOuterNormalIdx);
                    numPairOuterNormalIdx++;
                }
            }
        }
        else // intersection is on boundary
        {
            // fill the normal pair entries
            for(int pairIdx = 0; pairIdx < numPairs; ++pairIdx)
            {
                assert(pairData_[pairIdx].normalPair.second == -1);
                const auto distance = FacetCenter - elementCenter;
                this->pairData_[pairIdx].normalDistance = std::move(distance.two_norm());
                numPairOuterNormalIdx++;
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
        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(auto& data : pairData_)
        {
            // parallel Dofs and Distances
            data.parallelDofs.fill(-1);
            data.parallelDistances.fill(0.0);
        }

        // set basic global positions and stencil size definitions
        const auto& elementCenter = element_.geometry().center();

        // get the parallel Dofs
        const int parallelLocalIdx = intersection_.indexInInside();
        int numPairParallelIdx = 0;
        for(const auto& intersection : intersections(gridView_, element_))
        {
            if( facetIsNormal_(intersection.indexInInside(), parallelLocalIdx) )
            {
                // Store the parallel dimension of self cell in the direction of the axis
                auto parallelAxisIdx = directionIndex(intersection);
                pairData_[numPairParallelIdx].parallelDistances[0] = setParallelPairDistances_(element_, parallelAxisIdx);

                if( intersection.neighbor() )
                {
                    // If the normal intersection has a neighboring cell, go in and store the parallel information.
                    auto outerElement = intersection.outside();
                    pairData_[numPairParallelIdx].parallelDofs[0] = gridView_.indexSet().subIndex(outerElement, parallelLocalIdx, codimIntersection);
                    pairData_[numPairParallelIdx].parallelDistances[1] = setParallelPairDistances_(outerElement, parallelAxisIdx);
                }
                else  // No parallel neighbor available
                {
                    // If the intersection has no neighbor we have to deal with the virtual outer parallel dof
                    const auto& boundaryFacetCenter = intersection.geometry().center();
                    const auto distance = boundaryFacetCenter - elementCenter;
                    const auto virtualFirstParallelFaceDofPos = this->intersection_.geometry().center() + distance;

                    pairData_[numPairParallelIdx].virtualFirstParallelFaceDofPos = std::move(virtualFirstParallelFaceDofPos);
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
        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(auto& data : pairData_)
        {
            // parallel Dofs and Distances
            data.parallelDofs.fill(-1);
            data.parallelDistances.fill(0.0);
        }

        // set basic global positions and stencil size definitions
        unsigned int numParallelFaces = pairData_[0].parallelDistances.size();
        const auto& elementCenter = this->element_.geometry().center();

        // get the parallel Dofs
        const int parallelLocalIdx = intersection_.indexInInside();
        int numPairParallelIdx = 0;
        std::stack<Element> parallelElementStack;
        for(const auto& intersection : intersections(gridView_, element_))
        {
            if( facetIsNormal_(intersection.indexInInside(), parallelLocalIdx) )
            {
                if( intersection.neighbor() )
                {
                    auto parallelAxisIdx = directionIndex(intersection);
                    auto localNormalIntersectionIndex = intersection.indexInInside();
                    parallelElementStack.push(element_);

                    bool keepStacking =  (parallelElementStack.size() < numParallelFaces);
                    while(keepStacking)
                    {
                        auto e = parallelElementStack.top();
                        for(const auto& normalIntersection : intersections(gridView_, e))
                        {
                            if( normalIntersection.indexInInside() == localNormalIntersectionIndex )
                            {
                                if( normalIntersection.neighbor() )
                                {
                                    parallelElementStack.push(normalIntersection.outside());
                                    keepStacking = (parallelElementStack.size() < numParallelFaces);
                                }
                                else
                                {
                                keepStacking = false;
                                }
                            }
                        }
                    }

                    while(!parallelElementStack.empty())
                    {
                        if(parallelElementStack.size() > 1)
                        {
                            this->pairData_[numPairParallelIdx].parallelDofs[parallelElementStack.size()-2] = gridView_.indexSet().subIndex(parallelElementStack.top(), parallelLocalIdx, codimIntersection);
                        }
                        this->pairData_[numPairParallelIdx].parallelDistances[parallelElementStack.size()-1] = setParallelPairDistances_(parallelElementStack.top(), parallelAxisIdx);
                        parallelElementStack.pop();
                    }

                }
                else
                {
                    // If the intersection has no neighbor we have to deal with the virtual outer parallel dof
                    const auto& boundaryFacetCenter = intersection.geometry().center();

                    const auto distance = boundaryFacetCenter - elementCenter;
                    const auto virtualFirstParallelFaceDofPos = this->intersection_.geometry().center() + distance;

                    this->pairData_[numPairParallelIdx].virtualFirstParallelFaceDofPos = std::move(virtualFirstParallelFaceDofPos);

                    // The distance is saved doubled because with scvf.cellCenteredSelfToFirstParallelDistance
                    // an average between parallelDistances[0] and parallelDistances[1] will be computed
                    this->pairData_[numPairParallelIdx].parallelDistances[0] = std::move(distance.two_norm() * 2);
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

    //! Sets the information about the normal faces (within the element)
    void setNormalPairFirstInfo_(const int isIdx, const Element& element, const int numPairsIdx)
    {
        // store the inner normal dofIdx
        auto& dofIdx = pairData_[numPairsIdx].normalPair.first;
        dofIdx = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);

        // store the local normal facet index
        this->pairData_[numPairsIdx].localNormalFaceIdx = isIdx;
    }

    //! Sets the information about the normal faces (in the neighbor element)
    void setNormalPairSecondInfo_(const int isIdx, const Element& element, const int numPairsIdx)
    {
        // store the dofIdx
        auto& dofIdx = pairData_[numPairsIdx].normalPair.second;
        dofIdx = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);

        // store the element distance
        const auto& outerNormalFacet = getFacet_(isIdx, element);
        const auto outerNormalFacetPos = outerNormalFacet.geometry().center();
        const auto& innerNormalFacet = getFacet_(isIdx, element_);
        const auto innerNormalFacetPos = innerNormalFacet.geometry().center();
        this->pairData_[numPairsIdx].normalDistance = (innerNormalFacetPos - outerNormalFacetPos).two_norm();
    }

    //! Sets the information about the parallel distances
    Scalar setParallelPairDistances_(const Element& element, const int parallelAxisIdx) const
    {
        Scalar distance = 0;
        std::vector<GlobalPosition> faces;
        faces.reserve(numfacets);
        for (const auto& intElement : intersections(gridView_, element))
        {
            faces.push_back(intElement.geometry().center());
        }
        switch(parallelAxisIdx)
        {
            case 0:
            {
                distance = (faces[1] - faces[0]).two_norm();
            break;
            }
            case 1:
            {
                distance = (faces[3] - faces[2]).two_norm();
            break;
            }
            case 2:
            {
                assert(dim == 3);
                distance = (faces[5] - faces[4]).two_norm();
            break;
            }
            default:
            {
                DUNE_THROW(Dune::InvalidStateException, "Something went terribly wrong");
            }
        }
        return distance;
    }

    Intersection intersection_; //!< The intersection of interest
    const Element element_; //!< The respective element
    const typename Element::Geometry elementGeometry_; //!< Reference to the element geometry
    const GridView gridView_; //!< The grid view
    AxisData<Scalar, upwindSchemeOrder> axisData_; //!< Data related to forward and backward faces
    std::array<PairData<Scalar, GlobalPosition, upwindSchemeOrder>, numPairs> pairData_; //!< Collection of pair information related to normal and parallel faces
};

} // end namespace Dumux

#endif
