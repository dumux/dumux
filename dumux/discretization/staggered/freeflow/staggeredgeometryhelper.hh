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

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Data stored per sub face
 */
template<class Scalar, class GlobalPosition>
struct PairData
{
    int outerParallelFaceDofIdx;
    std::pair<int,int> normalPair;
    int localNormalFaceIdx;
    int globalCommonEntIdx;
    Scalar parallelDistance;
    Scalar normalDistance;
    GlobalPosition virtualOuterParallelFaceDofPos;
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
 *
 * A visualization of the variables that are in this document can be found in the following image:
 *
 * \image html staggeredgeometry.png
 *
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
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(intersection_.inside(), inIdx, codimIntersection);
    }

    /*!
     * \brief Returns the global dofIdx of the opposing intersection
     */
    int dofIndexOpposingFace() const
    {
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(this->intersection_.inside(), localOppositeIdx_(inIdx), codimIntersection);
    }

    /*!
     * \brief Returns the local index of the face (i.e. the intersection)
     */
    int localFaceIndex() const
    {
        return intersection_.indexInInside();
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
     * \brief Returns a copy of the pair data
     */
    auto pairData() const
    {
        return pairData_;
    }

     /*!
     * \brief Returns the dirction index of the facet (0 = x, 1 = y, 2 = z)
     */
    int directionIndex() const
    {
        return Dumux::directionIndex(std::move(intersection_.centerUnitOuterNormal()));
    }

private:
     /*!
     * \brief Fills all entries of the pair data
     */
    void fillPairData_()
    {
        const int localFacetIdx = intersection_.indexInInside();

        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(auto& data : pairData_)
        {
            data.outerParallelFaceDofIdx = -1;
            data.normalPair.second = -1;
            data.normalDistance = -1;
            data.parallelDistance = -1;
        }

        // set the inner parts of the normal pairs
        const auto localIndices = getLocalIndices_(localFacetIdx);
        setInnerIndices_(localIndices);

        // get the positions of the faces normal to the intersection within the element
        innerNormalFacePos_.reserve(numPairs);


        for(int i = 0; i < numPairs; ++i)
        {
           const auto& innerNormalFacet = getFacet_(localIndices.normalLocalDofIdx[i], element_);
           innerNormalFacePos_.emplace_back(innerNormalFacet.geometry().center());
        }

        // go into the direct neighbor element
        if(intersection_.neighbor())
        {
            // the direct neighbor element and the respective intersection index
            const auto& directNeighborElement = intersection_.outside();

            for(const auto& directNeighborElementIntersection : intersections(gridView_, directNeighborElement))
            {

                const int directNeighborElemIsIdx = directNeighborElementIntersection.indexInInside();
                // skip the directly neighboring face itself and its opposing one
                if(facetIsNormal_(directNeighborElemIsIdx, intersection_.indexInOutside()))
                {
                    setPairInfo_(directNeighborElemIsIdx, directNeighborElement, false);

                    // go into the adjacent neighbor element to get parallel dof info
                    if(directNeighborElementIntersection.neighbor())
                    {
                        const auto& diagonalNeighborElement = directNeighborElementIntersection.outside();
                        for(const auto& dIs : intersections(gridView_, diagonalNeighborElement))
                        {
                            if(facetIsNormal_(dIs.indexInInside(), directNeighborElementIntersection.indexInOutside()))
                                setPairInfo_(dIs.indexInInside(), diagonalNeighborElement, true);
                        }
                    }
                }
            }
        }
        else // intersection is on boundary
        {
            // find an intersection normal to the face
            for(const auto& normalIntersection : intersections(gridView_, element_))
            {
                if(facetIsNormal_(normalIntersection.indexInInside(), intersection_.indexInInside()) && normalIntersection.neighbor())
                {
                    const auto& neighborElement = normalIntersection.outside();

                    for(const auto& neighborElementIs : intersections(gridView_, neighborElement))
                    {
                        // iterate over facets sub-entities
                        if(facetIsNormal_(normalIntersection.indexInInside(), neighborElementIs.indexInInside()))
                            setPairInfo_(neighborElementIs.indexInInside(), neighborElement, true);
                    }
                }
            }

            // fill the normal pair entries
            for(int pairIdx = 0; pairIdx < numPairs; ++pairIdx)
            {
                assert(pairData_[pairIdx].normalPair.second == -1);
                const auto& elementCenter = this->element_.geometry().center();
                const auto& boundaryFacetCenter = intersection_.geometry().center();
                const auto distance = boundaryFacetCenter - elementCenter;

                pairData_[pairIdx].normalDistance = std::move(distance.two_norm());
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

    void setPairInfo_(const int isIdx, const Element& element, const bool isParallel)
    {
        const auto referenceElement = ReferenceElements::general(element_.geometry().type());

        // iterate over facets sub-entities
        for(int i = 0; i < numFacetSubEntities; ++i)
        {
            int localCommonEntIdx = referenceElement.subEntity(isIdx, 1, i, codimCommonEntity);
            int globalCommonEntIdx = localToGlobalCommonEntityIdx_(localCommonEntIdx, element);

            // fill the normal pair entries
            for(int pairIdx = 0; pairIdx < numPairs; ++pairIdx)
            {
                    if(globalCommonEntIdx == pairData_[pairIdx].globalCommonEntIdx)
                    {
                            auto& dofIdx = isParallel ? pairData_[pairIdx].outerParallelFaceDofIdx : pairData_[pairIdx].normalPair.second;
                            dofIdx = gridView_.indexSet().subIndex(element, isIdx, codimIntersection);
                            if(isParallel)
                            {
                                    const auto& selfFacet = getFacet_(intersection_.indexInInside(), element_);
                                    const auto& parallelFacet = getFacet_(isIdx, element);
                                    pairData_[pairIdx].parallelDistance = (selfFacet.geometry().center() - parallelFacet.geometry().center()).two_norm();
                            }
                            else
                            {
                                    const auto& outerNormalFacet = getFacet_(isIdx, element);
                                    const auto outerNormalFacetPos = outerNormalFacet.geometry().center();
                                    pairData_[pairIdx].normalDistance = (innerNormalFacePos_[pairIdx] - outerNormalFacetPos).two_norm();
                            }
                    }
            }
        }
    }

    template<class Indices>
    void setInnerIndices_(const Indices& indices)
    {
        for(int i = 0; i < numPairs; ++i)
        {
            this->pairData_[i].normalPair.first = this->gridView_.indexSet().subIndex(this->intersection_.inside(), indices.normalLocalDofIdx[i], codimIntersection);
            this->pairData_[i].globalCommonEntIdx = this->gridView_.indexSet().subIndex(this->intersection_.inside(), indices.localCommonEntIdx[i], codimCommonEntity);
            this->pairData_[i].localNormalFaceIdx = indices.normalLocalDofIdx[i];
        }
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
                const auto virtualOuterParallelFaceDofPos = this->intersection_.geometry().center() + distance;

                auto foundCorrectIdx = [otherIsIdx](const auto& x) { return x.localNormalFaceIdx == otherIsIdx; };
                const int index = std::find_if(this->pairData_.begin(), this->pairData_.end(), foundCorrectIdx) - this->pairData_.begin();

                this->pairData_[index].virtualOuterParallelFaceDofPos = std::move(virtualOuterParallelFaceDofPos);
                this->pairData_[index].parallelDistance = std::move(distance.two_norm());
            }
        }
    }

    auto getLocalIndices_(const int localFacetIdx) const
    {
        struct Indices
        {
            std::array<int, numPairs> normalLocalDofIdx;
            std::array<int, numPairs> localCommonEntIdx;
        };

        Indices indices;
        if (dim == 1)
            return indices;

        switch(localFacetIdx)
        {
            case 0:
                indices.normalLocalDofIdx[0] = 3;
                indices.normalLocalDofIdx[1] = 2;
                indices.localCommonEntIdx[0] = 2;
                indices.localCommonEntIdx[1] = 0;
                if(dim == 3)
                {
                    indices.normalLocalDofIdx[2] = 4;
                    indices.normalLocalDofIdx[3] = 5;
                    indices.localCommonEntIdx[2] = 4;
                    indices.localCommonEntIdx[3] = 8;
                }
                break;
            case 1:
                indices.normalLocalDofIdx[0] = 2;
                indices.normalLocalDofIdx[1] = 3;
                indices.localCommonEntIdx[0] = 1;
                indices.localCommonEntIdx[1] = 3;
                if(dim == 3)
                {
                    indices.normalLocalDofIdx[2] = 4;
                    indices.normalLocalDofIdx[3] = 5;
                    indices.localCommonEntIdx[2] = 5;
                    indices.localCommonEntIdx[3] = 9;
                }
                break;
            case 2:
                indices.normalLocalDofIdx[0] = 0;
                indices.normalLocalDofIdx[1] = 1;
                indices.localCommonEntIdx[0] = 0;
                indices.localCommonEntIdx[1] = 1;
                if(dim == 3)
                {
                    indices.normalLocalDofIdx[2] = 4;
                    indices.normalLocalDofIdx[3] = 5;
                    indices.localCommonEntIdx[2] = 6;
                    indices.localCommonEntIdx[3] = 10;
                }
                break;
            case 3:
                indices.normalLocalDofIdx[0] = 1;
                indices.normalLocalDofIdx[1] = 0;
                indices.localCommonEntIdx[0] = 3;
                indices.localCommonEntIdx[1] = 2;
                if(dim == 3)
                {
                    indices.normalLocalDofIdx[2] = 4;
                    indices.normalLocalDofIdx[3] = 5;
                    indices.localCommonEntIdx[2] = 7;
                    indices.localCommonEntIdx[3] = 11;
                }
                break;
            case 4:
                assert(dim == 3);
                indices.normalLocalDofIdx[0] = 0;
                indices.normalLocalDofIdx[1] = 1;
                indices.normalLocalDofIdx[2] = 3;
                indices.normalLocalDofIdx[3] = 2;
                indices.localCommonEntIdx[0] = 4;
                indices.localCommonEntIdx[1] = 5;
                indices.localCommonEntIdx[2] = 7;
                indices.localCommonEntIdx[3] = 6;
                break;
            case 5:
                assert(dim == 3);
                indices.normalLocalDofIdx[0] = 0;
                indices.normalLocalDofIdx[1] = 1;
                indices.normalLocalDofIdx[2] = 2;
                indices.normalLocalDofIdx[3] = 3;
                indices.localCommonEntIdx[0] = 8;
                indices.localCommonEntIdx[1] = 9;
                indices.localCommonEntIdx[2] = 10;
                indices.localCommonEntIdx[3] = 11;
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Something went terribly wrong");
        }
        return indices;
    }

    // TODO: check whether to use references here or not
    Intersection intersection_; //!< The intersection of interest
    const Element element_; //!< The respective element
    const typename Element::Geometry elementGeometry_; //!< Reference to the element geometry
    const GridView gridView_;
    std::array<PairData<Scalar, GlobalPosition>, numPairs> pairData_; //!< collection of pair information
    std::vector<GlobalPosition> innerNormalFacePos_;
};

} // end namespace Dumux

#endif
