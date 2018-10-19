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
 * \copydoc Dumux::FreeFlowStaggeredGeometryHelper
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GEOMETRY_HELPER_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/math.hh>
#include <type_traits>
#include <algorithm>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief Parallel Data stored per sub face
 */
template<class Scalar, class GlobalPosition>
struct PairData
{
    std::vector<int> parallelDofs;
    std::vector<Scalar> parallelDistances;
    std::pair<signed int,signed int> normalPair;
    int localNormalFaceIdx;
    Scalar normalDistance;
    GlobalPosition virtualOuterNormalFaceDofPos;
    GlobalPosition virtualFirstParallelFaceDofPos;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief In Axis Data stored per sub face
 */
template<class Scalar>
struct AxisData
{
    int selfDof;
    int oppositeDof;
    std::vector<int> inAxisForwardDofs;
    std::vector<int> inAxisBackwardDofs;
    Scalar selfToOppositeDistance;
    std::vector<Scalar> inAxisForwardDistances;
    std::vector<Scalar> inAxisBackwardDistances;
};

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
template<class GridView>
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
    static int order_;

public:

    FreeFlowStaggeredGeometryHelper(const Element& element, const GridView& gridView) : element_(element), elementGeometry_(element.geometry()), gridView_(gridView)
    { }

    //! update the local face
    template<class IntersectionMapper>
    void updateLocalFace(const IntersectionMapper& intersectionMapper_, const Intersection& intersection)
    {
        intersection_ = intersection;
        innerNormalFacePos_.clear();
        fillPairData_();
        fillAxisData_();
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

    //! \brief Returns the order needed by the scheme
    static int order()
    {
        return order_;
    }

    //! \brief Set the order needed by the scheme
    static void setOrder(const int order)
    {
        order_ = order;
    }

private:
    /*!
     * \brief Fills all entries of the in axis data
     */
    void fillAxisData_()
    {
        unsigned int numForwardBackwardAxisDofs = order_ - 1;
        const auto inIdx = intersection_.indexInInside();
        const auto oppIdx = localOppositeIdx_(inIdx);

        // Clear the containers before filling them
        this->axisData_.inAxisForwardDofs.clear();
        this->axisData_.inAxisBackwardDofs.clear();
        this->axisData_.inAxisForwardDistances.clear();
        this->axisData_.inAxisBackwardDistances.clear();

        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(int i = 0; i < numForwardBackwardAxisDofs; i++)
        {
            this->axisData_.inAxisForwardDofs.push_back(-1);
            this->axisData_.inAxisBackwardDofs.push_back(-1);
            this->axisData_.inAxisForwardDistances.push_back(0);
            this->axisData_.inAxisBackwardDistances.push_back(0);
        }

        // Set the self Dof
        this->axisData_.selfDof = gridView_.indexSet().subIndex(this->intersection_.inside(), inIdx, codimIntersection);
        const auto self = getFacet_(inIdx, element_);

        // Set the opposite Dof
        this->axisData_.oppositeDof = gridView_.indexSet().subIndex(this->intersection_.inside(), oppIdx, codimIntersection);
        const auto opposite = getFacet_(oppIdx, element_);

        // Set the Self to Opposite Distance
        this->axisData_.selfToOppositeDistance = (self.geometry().center() - opposite.geometry().center()).two_norm();

        // Set the Forward Dofs
        std::stack<Element> inAxisForwardElementStack;
        auto selfFace = getFacet_(inIdx, element_);
        if(intersection_.neighbor())
        {
            if (order_ > 1)
                inAxisForwardElementStack.push(intersection_.outside());
            bool keepStackingForward = (inAxisForwardElementStack.size() < order_-1);
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
                            keepStackingForward = (inAxisForwardElementStack.size() < order_-1);
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
            this->axisData_.inAxisForwardDofs[forwardIdx] = gridView_.indexSet().subIndex(inAxisForwardElementStack.top(), inIdx, codimIntersection);
            selfFace = getFacet_(inIdx, inAxisForwardElementStack.top());
            forwardFaceCoordinates[forwardIdx] = selfFace.geometry().center();
            inAxisForwardElementStack.pop();
        }
        for(int i = 0; i< forwardFaceCoordinates.size(); i++)
        {
            if(i == 0)
            {
                this->axisData_.inAxisForwardDistances[i] = (forwardFaceCoordinates[i] - self.geometry().center()).two_norm();
            }
            else
            {
                this->axisData_.inAxisForwardDistances[i] = (forwardFaceCoordinates[i]- forwardFaceCoordinates[i-1]).two_norm();
            }
        }

        // Set the Backward Dofs
        std::stack<Element> inAxisBackwardElementStack;
        auto oppFace = getFacet_(oppIdx, element_);
        for(const auto& intersection : intersections(gridView_, element_))
        {
            if(intersection.indexInInside() == oppIdx)
            {
                if(intersection.neighbor())
                {
                    if (order_ > 1)
                        inAxisBackwardElementStack.push(intersection.outside());
                    bool keepStackingBackward = (inAxisBackwardElementStack.size() < order_-1);
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
                                    keepStackingBackward = (inAxisBackwardElementStack.size() < order_-1);
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
     * \brief Fills all entries of the pair data
     */
    void fillPairData_()
    {
        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(auto& data : pairData_)
        {
            int numParallelDofs = order_;

            // parallel Dofs
            data.parallelDofs.clear();
            data.parallelDofs.resize(numParallelDofs, -1);

            // parallel Distances
            data.parallelDistances.clear();
            data.parallelDistances.resize(numParallelDofs + 1, 0.0);

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
                setNormalPairInfo_(innerElementIntersectionIdx, element_, numPairInnerNormalIdx, 0);
                numPairInnerNormalIdx++;

                innerNormalFacePos_.reserve(numPairs);
                const auto& innerNormalFacet = getFacet_(innerElementIntersection.indexInInside(), element_);
                innerNormalFacePos_.emplace_back(innerNormalFacet.geometry().center());
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
                    setNormalPairInfo_(directNeighborElemIsIdx, directNeighborElement, numPairOuterNormalIdx, 1);
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
                this->pairData_[pairIdx].virtualOuterNormalFaceDofPos = innerNormalFacePos_[pairIdx] + distance;
                this->pairData_[pairIdx].normalDistance = std::move(distance.two_norm());
                numPairOuterNormalIdx++;
            }
        }


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
                    bool keepStacking =  (parallelElementStack.size() < order_+1);
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
                                    keepStacking = (parallelElementStack.size() < order_+1);
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
                    const auto& elementCenter = this->element_.geometry().center();
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

    /*!
     * \brief Returns the global index of the common entity
     *
     * \param localIdx The local index of the common entity
     * \param element The element
     */
    int localToGlobalCommonEntityIdx_(const int localIdx, const Element& element) const
    {
        return this->gridView_.indexSet().subIndex(element, localIdx, codimCommonEntity);
    };

    auto getFacet_(const int localFacetIdx, const Element& element) const
    {
        return element.template subEntity <1> (localFacetIdx);
    };

    //! Sets the information about the normal faces
    void setNormalPairInfo_(const int isIdx, const Element& element, const int numPairsIdx, const int innerOuterIdx)
    {
        switch(innerOuterIdx)
        {
            case 0:
            {
                // store the inner normal dofIdx
                auto& dofIdx = pairData_[numPairsIdx].normalPair.first;
                dofIdx = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);

                // store the local normal facet index
                this->pairData_[numPairsIdx].localNormalFaceIdx = isIdx;
            break;
            }
            case 1:
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
            break;
            }
        }
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
    AxisData<Scalar> axisData_; //!< Data related to forward and backward faces
    std::array<PairData<Scalar, GlobalPosition>, numPairs> pairData_; //!< Collection of pair information related to normal and parallel faces
    std::vector<GlobalPosition> innerNormalFacePos_;
};

template<class GridView>
int FreeFlowStaggeredGeometryHelper<GridView>::order_ = 1;

} // end namespace Dumux

#endif
