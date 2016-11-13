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
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

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

        // initialize outerParallel (in boundary case, where there is no outer parallel dof)
        pairData_[0].outerParallel = 0;
        pairData_[1].outerParallel = 0;

        // set the inner parts of the normal pairs
        setInnerNormalPairs_(indexInInside);

        // go into the direct neighbor element
        if(intersection_.neighbor())
        {
            // the direct neighbor element and the respective intersection index
            const auto& directNeighbor = intersection_.outside();

            for(const auto& neighborIntersection : intersections(gridView_, directNeighbor))
            {
                // skip the directly neighboring face itself and its opposing one
                if(neighborIntersectionNormalSide_(neighborIntersection.indexInInside(), intersection_.indexInOutside()))
                {
                    // get facet
//                    const auto& facet = directNeighbor.template subEntity < 1 > (neighborIntersection.indexInInside());

                    // iterate over facets sub-entities // TODO: get number correctly
                    for(int i = 0; i < 2; ++i)
                    {
                        int localCommonEntIdx = referenceElement.subEntity(neighborIntersection.indexInInside(), 1, i, dim);
                        // const auto& commonEnt = directNeighbor.template subEntity < codimCommonEntity > (localCommonEntIdx);

//                         int globalCommonEntIdx = gridView_.indexSet().subIndex(directNeighbor, localCommonEntIdx, codimCommonEntity);
                        // int globalCommonEntIdx = localToGlobalEntityIdx(localCommonEntIdx);
                        int globalCommonEntIdx = localToGlobalEntityIdx_(localCommonEntIdx, directNeighbor);

//                         int dofIdx = gridView_.indexSet().index(commonEnt);

//                         if(globalCommonEntIdx != dofIdx)
//                             std::cout << "error!";

                        // std::cout << "pos: " << commonEnt.geometry().center() << std::endl;

                        if(globalCommonEntIdx == pairData_[0].globalCommonEntIdx)
                        {
                            pairData_[0].normalPair.second = gridView_.indexSet().subIndex(directNeighbor, neighborIntersection.indexInInside(), dim-1) + offset_;
                        }
                        if(globalCommonEntIdx == pairData_[1].globalCommonEntIdx)
                        {
                            pairData_[1].normalPair.second = gridView_.indexSet().subIndex(directNeighbor, neighborIntersection.indexInInside(), dim-1) + offset_;
                        }
                    }

                    // std::cout << "inters: " << intersection_.geometry().center() <<
                    // " || side : " << neighborIntersection.geometry().center() << std::endl;

                    // std::cout << "in element " << gridView_.indexSet().index(intersection_.inside()) << ", direct neighbor: " <<  gridView_.indexSet().index(directNeighbor)
                    //           << std::endl;


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

                                    // const auto& commonEnt = diagonalNeighbor.template subEntity < codimCommonEntity > (localCommonEntIdx);
                                    // std::cout << "globalCommEnt:" << commonEnt.geometry().center() << std::endl;
                                    //
                                    // int dofIdx = gridView_.indexSet().index(commonEnt);
                                    // if(globalCommonEntIdx != dofIdx)
                                    //                             std::cout << "error!";


                                    if(globalCommonEntIdx == pairData_[0].globalCommonEntIdx)
                                        pairData_[0].outerParallel = gridView_.indexSet().subIndex(diagonalNeighbor, dIs.indexInInside(), dim-1) + offset_;

                                    if(globalCommonEntIdx == pairData_[1].globalCommonEntIdx)
                                        pairData_[1].outerParallel = gridView_.indexSet().subIndex(diagonalNeighbor, dIs.indexInInside(), dim-1) + offset_;
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


    template<class G = GridView>
    typename std::enable_if<G::dimension == 2, void>::type
    setInnerNormalPairs_(const int directNeighborIsIdx)
    {
        int normalLocalDofIdx1;
        int normalLocalDofIdx2;
        int localCommonEntIdx1;
        int localCommonEntIdx2;

        switch(directNeighborIsIdx)
        {
            case 0:
                normalLocalDofIdx1 = 3;
                normalLocalDofIdx2 = 2;
                localCommonEntIdx1 = 2;
                localCommonEntIdx2 = 0;
                break;
            case 1:
                normalLocalDofIdx1 = 2;
                normalLocalDofIdx2 = 3;
                localCommonEntIdx1 = 1;
                localCommonEntIdx2 = 3;
                break;
            case 2:
                normalLocalDofIdx1 = 0;
                normalLocalDofIdx2 = 1;
                localCommonEntIdx1 = 0;
                localCommonEntIdx2 = 1;
                break;
            case 3:
                normalLocalDofIdx1 = 1;
                normalLocalDofIdx2 = 0;
                localCommonEntIdx1 = 3;
                localCommonEntIdx2 = 2;
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "Something went terribly wrong");
        }
        pairData_[0].normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), normalLocalDofIdx1, dim-1) + offset_;
        pairData_[1].normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), normalLocalDofIdx2, dim-1) + offset_;
        pairData_[0].globalCommonEntIdx = gridView_.indexSet().subIndex(intersection_.inside(), localCommonEntIdx1, codimCommonEntity);
        pairData_[1].globalCommonEntIdx = gridView_.indexSet().subIndex(intersection_.inside(), localCommonEntIdx2, codimCommonEntity);
    }

    template<class G = GridView>
    typename std::enable_if<G::dimension == 3, void>::type
    setInnerNormalPairs_(const int directNeighborIsIdx)
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
