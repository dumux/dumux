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

namespace Dumux
{

template<class Scalar>
struct PairData
{
    int outerParallel;
    std::pair<int,int> normalPair;
    int globalCommonEntIdx;
    Scalar parallelDistance;
    Scalar normalDistance;
};


//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
// class StaggeredGeometryHelper<GridView, dim>
class StaggeredGeometryHelper
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


public:


    StaggeredGeometryHelper(const Intersection& intersection, const GridView& gridView)
    : intersection_(intersection), element_(intersection.inside()), elementGeometry_(element_.geometry()), gridView_(gridView), offset_(gridView.size(0))
    {
        fillPairData();
    }

    int dofIdxSelf() const
    {
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(intersection_.inside(), inIdx, dim-1) + offset_;
    }

    int dofIdxOpposite() const
    {
        //TODO: use proper intersection mapper!
        const auto inIdx = intersection_.indexInInside();
        return gridView_.indexSet().subIndex(intersection_.inside(), localOppositeIdx_(inIdx), dim-1) + offset_;
    }


    auto pairData() const
    {
        return pairData_;
    }


    void fillPairData()
    {
        const auto& referenceElement = ReferenceElements::general(element_.geometry().type());
        const int indexInInside = intersection_.indexInInside();

        // initialize values that could remain unitialized if the intersection lies on a boundary
        for(auto& data : pairData_)
        {
            data.outerParallel = 0;
            data.normalDistance = 0;
            data.parallelDistance = 0;
        }


        // set the inner parts of the normal pairs
        const auto localInnerNormalDofIndices = getLocalInnerNormalDofIndices_(indexInInside);
        setInnerNormalPairs_(localInnerNormalDofIndices);

        // get the positions of the faces normal to the intersection within the element
        std::vector<GlobalPosition> innerNormalFacePos;
        innerNormalFacePos.reserve(numPairs);
        if(dimWorld == 2)
        {
            const auto& innerNormalFacet1 = element_.template subEntity <1> (localInnerNormalDofIndices.normalLocalDofIdx1);
            const auto& innerNormalFacet2 = element_.template subEntity <1> (localInnerNormalDofIndices.normalLocalDofIdx2);
            innerNormalFacePos.emplace_back(innerNormalFacet1.geometry().center());
            innerNormalFacePos.emplace_back(innerNormalFacet2.geometry().center());
        }
        if(dimWorld == 3)
        {
            DUNE_THROW(Dune::NotImplemented, "3d not ready yet");
        }


        // go into the direct neighbor element
        if(intersection_.neighbor())
        {
            // the direct neighbor element and the respective intersection index
            const auto& directNeighbor = intersection_.outside();

            for(const auto& neighborIntersection : intersections(gridView_, directNeighbor))
            {
                const int neighborIsIdx = neighborIntersection.indexInInside();
                // skip the directly neighboring face itself and its opposing one
                if(neighborIntersectionNormalSide_(neighborIsIdx, intersection_.indexInOutside()))
                {
                    // iterate over facets sub-entities // TODO: get number correctly
                    for(int i = 0; i < 2; ++i)
                    {
                        int localCommonEntIdx = referenceElement.subEntity(neighborIsIdx, 1, i, dim);
                        int globalCommonEntIdx = localToGlobalEntityIdx_(localCommonEntIdx, directNeighbor);

                        // fill the normal pair entries
                        for(int pairIdx = 0; pairIdx < numPairs; ++pairIdx)
                        {
                            if(globalCommonEntIdx == pairData_[pairIdx].globalCommonEntIdx)
                            {
                                pairData_[pairIdx].normalPair.second = gridView_.indexSet().subIndex(directNeighbor, neighborIsIdx, dim-1) + offset_;
                                const auto& outerNormalFacet = directNeighbor.template subEntity <1> (neighborIsIdx);
                                const auto outerNormalFacetPos = outerNormalFacet.geometry().center();
                                pairData_[pairIdx].normalDistance = (innerNormalFacePos[pairIdx] - outerNormalFacetPos).two_norm();
                            }
                        }
                    }



                    // go into the adjacent neighbor element
                    if(neighborIntersection.neighbor())
                    {
                        const auto& diagonalNeighbor = neighborIntersection.outside();
                        for(const auto& dIs : intersections(gridView_, diagonalNeighbor))
                        {
                            if(neighborIntersectionNormalSide_(dIs.indexInInside(), neighborIntersection.indexInOutside()))
                            {
                                for(int i = 0; i < 2; ++i)
                                {
                                    int localCommonEntIdx = referenceElement.subEntity(dIs.indexInInside(), 1, i, dim);
                                    int globalCommonEntIdx = localToGlobalEntityIdx_(localCommonEntIdx, diagonalNeighbor);


                                    for(int pairIdx = 0; pairIdx < numPairs; ++pairIdx)
                                    {
                                        if(globalCommonEntIdx == pairData_[pairIdx].globalCommonEntIdx)
                                        {
                                            pairData_[pairIdx].outerParallel = gridView_.indexSet().subIndex(diagonalNeighbor, dIs.indexInInside(), dim-1) + offset_;
                                            const auto& selfFacet = element_.template subEntity <1> (indexInInside);
                                            const auto& parallelFacet = diagonalNeighbor.template subEntity <1> (dIs.indexInInside());
                                            pairData_[pairIdx].parallelDistance = (selfFacet.geometry().center() - parallelFacet.geometry().center()).two_norm();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

private:

    int localOppositeIdx_(const int idx) const
    {
        return (idx % 2) ? (idx - 1) : (idx + 1);
    }

    bool neighborIntersectionNormalSide_(const int isIdx, const int neighborIsIdx) const
    {
        return !(isIdx == neighborIsIdx || localOppositeIdx_(isIdx) == neighborIsIdx);
    };

    int localToGlobalEntityIdx_(const int localIdx, const Element& element)
    {
        return gridView_.indexSet().subIndex(element, localIdx, codimCommonEntity);
    };




    // specializations for 2D ***************************************************************************************************
    template<class G = GridView, typename std::enable_if<G::dimension == 2, int>::type = 0>
    auto getLocalInnerNormalDofIndices_(const int directNeighborIsIdx)
    {
        struct Indices
        {
            int normalLocalDofIdx1;
            int normalLocalDofIdx2;
            int localCommonEntIdx1;
            int localCommonEntIdx2;
        };

        Indices indices;

        switch(directNeighborIsIdx)
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

    template<class Indices, class G = GridView>
    typename std::enable_if<G::dimension == 2, void>::type
    setInnerNormalPairs_(const Indices& indices)
    {
        pairData_[0].normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), indices.normalLocalDofIdx1, dim-1) + offset_;
        pairData_[1].normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), indices.normalLocalDofIdx2, dim-1) + offset_;
        pairData_[0].globalCommonEntIdx = gridView_.indexSet().subIndex(intersection_.inside(), indices.localCommonEntIdx1, codimCommonEntity);
        pairData_[1].globalCommonEntIdx = gridView_.indexSet().subIndex(intersection_.inside(), indices.localCommonEntIdx2, codimCommonEntity);
    }


    template<class G = GridView, typename std::enable_if<G::dimension == 3, int>::type = 0>
    auto getLocalInnerNormalDofIndices_(const int directNeighborIsIdx)
    {
        struct Indices
        {
            int normalLocalDofIdx1;
            int normalLocalDofIdx2;
            int normalLocalDofIdx3;
            int normalLocalDofIdx4;
            int localCommonEntIdx1;
            int localCommonEntIdx2;
            int localCommonEntIdx3;
            int localCommonEntIdx4;
        };

        Indices indices;

        switch(directNeighborIsIdx)
        {
            default:
                DUNE_THROW(Dune::NotImplemented, "3d helper not ready yet");
        }
        return indices;
    }

    template<class Indices, class G = GridView>
    typename std::enable_if<G::dimension == 3, void>::type
    setInnerNormalPairs_(const Indices& indices)
    {
        // TODO: 3D
        DUNE_THROW(Dune::NotImplemented, "3d helper not ready yet");
    }


    const Intersection& intersection_;
    const Element& element_;
    const typename Element::Geometry& elementGeometry_; //! Reference to the element geometry
    const GridView gridView_;
    const int offset_;

    std::array<PairData<Scalar>, numPairs> pairData_;


};



} // end namespace Dumux

#endif
