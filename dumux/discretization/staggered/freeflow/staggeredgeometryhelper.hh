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
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the staggered discretization method
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

template<class Scalar, class GlobalPosition>
struct PairData
{
    int outerParallelFaceDofIdx;
    std::pair<int,int> normalPair;
    int localNormalFaceIdx;
    int globalCommonEntIdx;
    Scalar parallelDistance;
    Scalar normalDistance;
    GlobalPosition virtualOuterNormalFaceDofPos;
    GlobalPosition virtualOuterParallelFaceDofPos;
};


 /*!
 * \brief Returns the dirction index of the facet (0 = x, 1 = y, 2 = z)
 */
template<class Vector>
inline static int directionIndex(Vector&& vector)
{
    const auto eps = 1e-8;
    const int idx = std::find_if(vector.begin(), vector.end(), [eps](const auto& x) { return std::abs(x) > eps; } ) - vector.begin();
    return idx;
}

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, int dim = GridView::dimension >
class StaggeredGeometryHelper
{};

template<class GridView>
class BaseStaggeredGeometryHelper
{
    using Scalar = typename GridView::ctype;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr int numPairs = (dimWorld == 2) ? 2 : 4;

    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using PointVector = std::vector<GlobalPosition>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

    //TODO include assert that checks for quad geometry
    static constexpr int codimCommonEntity = 2; //TODO: 3d?
    static constexpr int numFacetSubEntities = 2; // TODO: 3d?

    using Implementation = typename Dumux::StaggeredGeometryHelper<GridView, dim>;

public:
    BaseStaggeredGeometryHelper(const Intersection& intersection, const GridView& gridView)
    : intersection_(intersection), element_(intersection.inside()), elementGeometry_(element_.geometry()), gridView_(gridView)
    {
        fillPairData_();
    }

     /*!
     * \brief Returns the global dofIdx of the intersection itself
     */
    int dofIdxSelf() const
    {
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(intersection_.inside(), inIdx, dim-1);
    }

     /*!
     * \brief Returns the global dofIdx of the opposing intersection
     */
    int dofIdxOpposite() const
    {
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(this->intersection_.inside(), localOppositeIdx_(inIdx), dim-1);
    }

     /*!
     * \brief Returns the distance between dofSelf and dofOpposite
     */
    Scalar selfToOppositeDistance() const
    {
        const auto inIdx = intersection_.indexInInside();
        const auto self = element_.template subEntity <1> (inIdx);
        const auto opposite = element_.template subEntity <1> (localOppositeIdx_(inIdx));
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
        const auto localIndices = asImp_().getLocalIndices_(localFacetIdx);
        asImp_().setInnerIndices_(localIndices);

        // get the positions of the faces normal to the intersection within the element
        innerNormalFacePos_.reserve(numPairs);
        if(dimWorld == 2)
        {
            const auto& innerNormalFacet1 = getFacet_(localIndices.normalLocalDofIdx1, element_);
            const auto& innerNormalFacet2 = getFacet_(localIndices.normalLocalDofIdx2, element_);
            innerNormalFacePos_.emplace_back(innerNormalFacet1.geometry().center());
            innerNormalFacePos_.emplace_back(innerNormalFacet2.geometry().center());
        }
//        if(dimWorld == 3)
//        {
////            DUNE_THROW(Dune::NotImplemented, "3d not ready yet");
//            const auto& innerNormalFacet1 = getFacet_(localIndices.normalLocalDofIdx1, element_);
//            const auto& innerNormalFacet2 = getFacet_(localIndices.normalLocalDofIdx2, element_);
//            const auto& innerNormalFacet3 = getFacet_(localIndices.normalLocalDofIdx3, element_);
//            const auto& innerNormalFacet4 = getFacet_(localIndices.normalLocalDofIdx4, element_);
//            innerNormalFacePos_.emplace_back(innerNormalFacet1.geometry().center());
//            innerNormalFacePos_.emplace_back(innerNormalFacet2.geometry().center());
//            innerNormalFacePos_.emplace_back(innerNormalFacet3.geometry().center());
//            innerNormalFacePos_.emplace_back(innerNormalFacet4.geometry().center());
//        }

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

                pairData_[pairIdx].virtualOuterNormalFaceDofPos = innerNormalFacePos_[pairIdx] + distance;
                pairData_[pairIdx].normalDistance = std::move(distance.two_norm());
            }
        }

        asImp_().treatVirtualOuterParallelFaceDofs_();
    }

protected:
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
        const int numFacetSubEntities = 2;
        const auto& referenceElement = ReferenceElements::general(element_.geometry().type());

        // iterate over facets sub-entities
        for(int i = 0; i < numFacetSubEntities; ++i)
        {
            int localCommonEntIdx = referenceElement.subEntity(isIdx, 1, i, dim);
            int globalCommonEntIdx = localToGlobalCommonEntityIdx_(localCommonEntIdx, element);

            // fill the normal pair entries
            for(int pairIdx = 0; pairIdx < numPairs; ++pairIdx)
            {
                    if(globalCommonEntIdx == pairData_[pairIdx].globalCommonEntIdx)
                    {
                            auto& dofIdx = isParallel ? pairData_[pairIdx].outerParallelFaceDofIdx : pairData_[pairIdx].normalPair.second;
                            dofIdx = gridView_.indexSet().subIndex(element, isIdx, dim-1);
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

    // TODO: check whether to use references here or not
    const Intersection intersection_; //! The intersection of interest
    const Element element_; //! The respective element
    const typename Element::Geometry elementGeometry_; //! Reference to the element geometry
    const GridView gridView_;
    std::array<PairData<Scalar, GlobalPosition>, numPairs> pairData_; //! collection of pair information

    std::vector<GlobalPosition> innerNormalFacePos_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};


template<class GridView>
class StaggeredGeometryHelper<GridView, 2> : public BaseStaggeredGeometryHelper<GridView>
{
    friend class BaseStaggeredGeometryHelper<GridView>;
    using Scalar = typename GridView::ctype;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr int numPairs = (dimWorld == 2) ? 2 : 4;

    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using PointVector = std::vector<GlobalPosition>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

    //TODO include assert that checks for quad geometry
    static constexpr int codimCommonEntity = 2; //TODO: 3d?

    using ParentType = BaseStaggeredGeometryHelper<GridView>;

public:
    StaggeredGeometryHelper(const Intersection& intersection, const GridView& gridView)
    : ParentType(intersection, gridView)
    {}

private:
    static auto getLocalIndices_(const int localFacetIdx)
    {
        struct Indices
        {
            int normalLocalDofIdx1;
            int normalLocalDofIdx2;
            int localCommonEntIdx1;
            int localCommonEntIdx2;
        };

        Indices indices;

        switch(localFacetIdx)
        {
            case 0:
                indices.normalLocalDofIdx1 = 3;
                indices.normalLocalDofIdx2 = 2;
                indices.localCommonEntIdx1 = 2;
                indices.localCommonEntIdx2 = 0;
                break;
            case 1:
                indices.normalLocalDofIdx1 = 2;
                indices.normalLocalDofIdx2 = 3;
                indices.localCommonEntIdx1 = 1;
                indices.localCommonEntIdx2 = 3;
                break;
            case 2:
                indices.normalLocalDofIdx1 = 0;
                indices.normalLocalDofIdx2 = 1;
                indices.localCommonEntIdx1 = 0;
                indices.localCommonEntIdx2 = 1;
                break;
            case 3:
                indices.normalLocalDofIdx1 = 1;
                indices.normalLocalDofIdx2 = 0;
                indices.localCommonEntIdx1 = 3;
                indices.localCommonEntIdx2 = 2;
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Something went terribly wrong");
        }
        return indices;
    }

    template<class Indices>
    void setInnerIndices_(const Indices& indices)
    {
        this->pairData_[0].normalPair.first = this->gridView_.indexSet().subIndex(this->intersection_.inside(), indices.normalLocalDofIdx1, dim-1);
        this->pairData_[1].normalPair.first = this->gridView_.indexSet().subIndex(this->intersection_.inside(), indices.normalLocalDofIdx2, dim-1);
        this->pairData_[0].globalCommonEntIdx = this->gridView_.indexSet().subIndex(this->intersection_.inside(), indices.localCommonEntIdx1, codimCommonEntity);
        this->pairData_[1].globalCommonEntIdx = this->gridView_.indexSet().subIndex(this->intersection_.inside(), indices.localCommonEntIdx2, codimCommonEntity);

        this->pairData_[0].localNormalFaceIdx = indices.normalLocalDofIdx1;
        this->pairData_[1].localNormalFaceIdx = indices.normalLocalDofIdx2;
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
            // check if any of these intersection, which are normal to our facet, lie on the domain boundary
            if(( !is.neighbor() ) &&  this->facetIsNormal_(localIntersectionIdx, is.indexInInside()) )
            {
                const auto& elementCenter = this->element_.geometry().center();
                const auto& boundaryFacetCenter = is.geometry().center();

                const auto distance = boundaryFacetCenter - elementCenter;
                const auto virtualOuterParallelFaceDofPos = this->intersection_.geometry().center() + distance;

                int index;
                switch(localIntersectionIdx)
                {
                    case 0:
                        index = (is.indexInInside() == 3) ? 0 : 1;
                        break;
                    case 1:
                        index = (is.indexInInside() == 2) ? 0 : 1;
                        break;
                    case 2:
                        index = (is.indexInInside() == 0) ? 0 : 1;
                        break;
                    case 3:
                        index = (is.indexInInside() == 1) ? 0 : 1;
                        break;
                    default:
                        DUNE_THROW(Dune::InvalidStateException, "Something went terribly wrong");
                }

                assert(this->pairData_[index].outerParallelFaceDofIdx == -1);

                this->pairData_[index].virtualOuterParallelFaceDofPos = std::move(virtualOuterParallelFaceDofPos);
                this->pairData_[index].parallelDistance = std::move(distance.two_norm());
            }
        }
    }
};

template<class GridView>
class StaggeredGeometryHelper<GridView, 3> : public BaseStaggeredGeometryHelper<GridView>
{
    friend class BaseStaggeredGeometryHelper<GridView>;
    using Scalar = typename GridView::ctype;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using PointVector = std::vector<GlobalPosition>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;


    //TODO include assert that checks for quad geometry
    static constexpr int codimCommonEntity = 2; //TODO: 3d?

    using ParentType = BaseStaggeredGeometryHelper<GridView>;

public:
    StaggeredGeometryHelper(const Intersection& intersection, const GridView& gridView)
    : ParentType(intersection, gridView)
    {}

private:
    auto getLocalIndices_(const int directdirectNeighborElemIsIdx)
    {
        // TODO: 3D
        DUNE_THROW(Dune::NotImplemented, "3d helper not ready yet");
    }

    template<class Indices>
    void setInnerIndices_(const Indices& indices)
    {
        // TODO: 3D
        DUNE_THROW(Dune::NotImplemented, "3d helper not ready yet");
    }
};



} // end namespace Dumux

#endif
