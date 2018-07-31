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
 * \brief Data stored per sub face
 */
template<class Scalar, class GlobalPosition>
struct PairData
{
    int firstParallelFaceDofIdx;
    int secondParallelFaceDofIdx;
    std::pair<int,int> normalPair;
    int localNormalFaceIdx;
    Scalar selfParallelElementDistance;
    Scalar firstParallelElementDistance;
    Scalar secondParallelElementDistance;
    Scalar normalDistance;
    GlobalPosition virtualOuterNormalFaceDofPos;
    GlobalPosition virtualFirstParallelFaceDofPos;

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
    }

    /*!
     * \brief Returns the global dofIdx of the intersection itself
     */
    int dofIndex() const
    {
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(this->intersection_.inside(), inIdx, codimIntersection);
    }

    /*!
     * \brief Returns the global dofIdx of the opposing intersection
     */
    int dofIndexOpposingFace() const
    {
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(this->intersection_.inside(), localOppositeIdx_(inIdx), codimIntersection);
    }

    /*!
     * \brief Returns the global dofIdx of the next intersection
     */
    int dofIndexPreviousFace() const
    {
        const int inIdx = intersection_.indexInInside();
        const int oppIdx = localOppositeIdx_(inIdx);
        int previousIdx = -1;
        for(const auto& intersection1 : intersections(gridView_, element_))
        {
            if( (intersection1.indexInInside() == oppIdx ) )
            {
                if( intersection1.neighbor() )
                {
                    const auto directPreviousElement = intersection1.outside();
                    previousIdx = gridView_.indexSet().subIndex(directPreviousElement, oppIdx, codimIntersection);
                }
            }
        }
        return previousIdx;
    }

    /*!
     * \brief Returns the global dofIdx of the next intersection
     */
    int dofIndexForwardFace() const
    {
        int dofIndexForwardFace = -1;
        if(intersection_.neighbor())
        {
            // the forward neighbor element and the respective intersection index
            const auto& directNeighborElement = intersection_.outside();              // Stores the forward neighboring element
            const int directForwardElemIsIdx = intersection_.indexInOutside();        // Stores the current scvf index on the outside neighbor cell
            // pass the forward element, and the opposite of the original facet's index in the forward cell to the subindex function.
            dofIndexForwardFace = gridView_.indexSet().subIndex(directNeighborElement, localOppositeIdx_(directForwardElemIsIdx), codimIntersection);
        }
        return dofIndexForwardFace;
    }

    /*!
     * \brief Returns the local index of the face (i.e. the intersection)
     */
    int localFaceIndex() const
    {
        return intersection_.indexInInside();
    }

    /*!
     * \brief Returns the distance between dofOpposite dofPrevious
     */
    Scalar oppositeToPreviousDistance() const
    {
        const auto inIdx = intersection_.indexInInside();
        Scalar oppositeToPreviousDistance = 0.0;
        for(const auto& intersectionSelf : intersections(gridView_, element_))
        {
            if( intersectionSelf.indexInInside() == localOppositeIdx_(inIdx) )
            {
                if( intersectionSelf.neighbor() )
                {
                    const auto oppOutsideIdx = intersectionSelf.indexInOutside();
                    const auto& previousEnt = intersectionSelf.outside();
                    const auto previousIdx = localOppositeIdx_(oppOutsideIdx);
                    const auto previousFacet = getFacet_(previousIdx, previousEnt);
                    oppositeToPreviousDistance = (intersectionSelf.geometry().center() - previousFacet.geometry().center()).two_norm();
                }
            }
        }
        return oppositeToPreviousDistance;
    }

    /*!
     * \brief Returns the distance between dofSelf and dofOpposite
     */
    Scalar selfToOppositeDistance() const
    {
        const auto inIdx = intersection_.indexInInside();
        const auto self = getFacet_(inIdx, element_);
        const auto opposite = getFacet_(localOppositeIdx_(inIdx), element_);
        return (self.geometry().center() - opposite.geometry().center()).two_norm();
    }

    /*!
     * \brief Returns the distance between dofSelf and dofOpposite
     */
    Scalar forwardToSelfDistance() const
    {
        const auto inIdx = intersection_.indexInInside();
        Scalar forwardToSelfDistance = 0.0;
        if(intersection_.neighbor())
        {
          const auto& selfFacet = getFacet_(inIdx, element_);
          const auto selfOutIdx = intersection_.indexInOutside();
          const auto forwardIdx = localOppositeIdx_(selfOutIdx);
          const auto& forwardEnt = intersection_.outside();
          const auto forwardFacet = getFacet_(forwardIdx, forwardEnt);
          forwardToSelfDistance = (forwardFacet.geometry().center() - selfFacet.geometry().center()).two_norm();
        }
        return forwardToSelfDistance;
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
     * \brief Fills all entries of the pair data
     */
    void fillPairData_()
    {
        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(auto& data : pairData_)
        {
            data.normalPair.second = -1;
            data.firstParallelFaceDofIdx = -1;
            data.secondParallelFaceDofIdx = -1;
            data.selfParallelElementDistance = 0.0;
            data.firstParallelElementDistance = 0.0;
            data.secondParallelElementDistance = 0.0;
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
            }
        }

        // get the parallel Dofs
        const int parallelIdx = intersection_.indexInInside();
        int numPairParallelIdx = 0;
        for(const auto& Intersection1 : intersections(gridView_, element_))
        {
            if( facetIsNormal_(Intersection1.indexInInside(), parallelIdx) )
            {
                if( Intersection1.neighbor() )
                {
                    auto parallelAxisIdx = directionIndex(Intersection1);
                    setParallelPairInfo_(intersection_.indexInInside(), element_, 0, numPairParallelIdx, parallelAxisIdx);
                    const auto& parallelElement1 = Intersection1.outside();
                    setParallelPairInfo_(parallelIdx, parallelElement1, 1, numPairParallelIdx, parallelAxisIdx);
                    for(const auto& Intersection2 : intersections(gridView_, parallelElement1) )
                    {
                        if( facetIsNormal_(Intersection2.indexInInside(), parallelIdx)
                            && ( (Intersection2.geometry().center() - elementCenter).two_norm() > (Intersection1.geometry().center() - elementCenter).two_norm() ) )
                        {
                            if( Intersection2.neighbor() )
                            {
                                const auto& parallelElement2 = Intersection2.outside();
                                setParallelPairInfo_(parallelIdx, parallelElement2, 2, numPairParallelIdx, parallelAxisIdx);
                            }
                        }
                    }
                }
                numPairParallelIdx++;
            }
        }
        treatVirtualOuterParallelFaceDofs_();
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

    void setParallelPairInfo_(const int isIdx, const Element& element, const int parallelDegreeIdx, const int numPairsIdx, const int parallelAxisIdx)
    {
        switch(parallelDegreeIdx)
        {
            case 0:
            {
                  //No need to store the dofIdx, is already stored as dofIdx(Self)
                  // store the self element distance
                  std::vector<GlobalPosition> facetSelf;
                  facetSelf.reserve(numfacets);
                  for(const auto& Int1 : intersections(gridView_, element_))
                  {
                      facetSelf.push_back(Int1.geometry().center());
                  }
                  this->pairData_[numPairsIdx].firstParallelElementDistance = setParallelPairDistances(facetSelf, parallelAxisIdx);
            break;
            }
            case 1:
            {
                  // store the dofIdx
                  auto& dofIdx = pairData_[numPairsIdx].firstParallelFaceDofIdx;
                  dofIdx = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);

                  // store the first element distance
                  std::vector<GlobalPosition> facetFirst;
                  facetFirst.reserve(numfacets);
                  for(const auto& Int1 : intersections(gridView_, element))
                  {
                      facetFirst.push_back(Int1.geometry().center());
                  }
                  this->pairData_[numPairsIdx].firstParallelElementDistance = setParallelPairDistances(facetFirst, parallelAxisIdx);
            break;
            }
            case 2:
            {
                  // store the dofIdx
                  auto& dofIdx = pairData_[numPairsIdx].secondParallelFaceDofIdx;
                  dofIdx = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);

                  // store the element distance
                  std::vector<GlobalPosition> facetSecond;
                  facetSecond.reserve(numfacets);
                  for(const auto& Int2 : intersections(gridView_, element))
                  {
                      facetSecond.push_back(Int2.geometry().center());
                  }
                  this->pairData_[numPairsIdx].secondParallelElementDistance = setParallelPairDistances(facetSecond, parallelAxisIdx);
            break;
            }
            default:
            {
                DUNE_THROW(Dune::InvalidStateException, "Something went terribly wrong");
            }
        }
    }

    Scalar setParallelPairDistances(const std::vector<GlobalPosition> faces, const int parallelAxisIdx)
    {
        Scalar distance = 0;
        if(parallelAxisIdx == 0)
        {
            distance = (faces[1] - faces[0]).two_norm();
        }
        else if(parallelAxisIdx == 1)
        {
            distance = (faces[3] - faces[2]).two_norm();
        }
        else
        {
            static_assert(dim==3, "3D face indices are called for a non-3D simulation")
            distance = (faces[5] - faces[4]).two_norm();
        }
        return distance;
    }

   /*!
    * \brief Sets the position of "virtual" dofs on facets which are perpendicular to facets on the domain boundary.
    *        This is required to account e.g. for wall friction.
    */
    void treatVirtualOuterParallelFaceDofs_()
    {
        // get the local index of the facet we are dealing with in this class
        const int localIntersectionIdx = this->intersection_.indexInInside();

        // iterate over all intersections of the element, this facet is part of
        for(const auto& is : intersections(this->gridView_, this->element_))
        {
            const int otherIsIdx = is.indexInInside();
            // check if any of these intersection, which are normal to our facet, lie on the domain boundary
            if(( !is.neighbor() ) &&  this->facetIsNormal_(localIntersectionIdx, otherIsIdx) )
            {
                const auto& elementCenter = this->element_.geometry().center();
                const auto& boundaryFacetCenter = is.geometry().center();

                const auto distance = boundaryFacetCenter - elementCenter;
                const auto virtualFirstParallelFaceDofPos = this->intersection_.geometry().center() + distance;

                auto foundCorrectIdx = [otherIsIdx](const auto& x) { return x.localNormalFaceIdx == otherIsIdx; };
                const int index = std::find_if(this->pairData_.begin(), this->pairData_.end(), foundCorrectIdx) - this->pairData_.begin();

                this->pairData_[index].virtualFirstParallelFaceDofPos = std::move(virtualFirstParallelFaceDofPos);
                this->pairData_[index].firstParallelElementDistance = std::move(distance.two_norm());
            }
        }
    }

    Intersection intersection_; //!< The intersection of interest
    const Element element_; //!< The respective element
    const typename Element::Geometry elementGeometry_; //!< Reference to the element geometry
    const GridView gridView_;
    std::array<PairData<Scalar, GlobalPosition>, numPairs> pairData_; //!< collection of pair information
    std::vector<GlobalPosition> innerNormalFacePos_;
};

} // end namespace Dumux

#endif
