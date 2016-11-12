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
 *        for the staggered discretizazion method
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
    std::pair<int,int> parallelPair;
    std::pair<int,int> normalPair;
    int localCommonEntIdx;
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


public:


    StaggeredGeometryHelper(const Intersection& intersection, const GridView& gridView)
    : intersection_(intersection), elementGeometry_(intersection.inside().geometry()), gridView_(gridView), offset_(gridView.size(0))
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
        const auto element = intersection_.inside();
        const auto& referenceElement = ReferenceElements::general(element.geometry().type());

        const int indexInInside = intersection_.indexInInside();


        // store the element's own face info first
        pairData_[0].parallelPair.first = dofIdxSelf();
        pairData_[1].parallelPair.first = dofIdxSelf();

        setInnerNormalPairs_(indexInInside);



        // go into the direct neigbor element
        if(intersection_.neighbor())
        {
            const int neighborIsIdx = intersection_.indexInOutside();

            const auto& directNeighbor = intersection_.outside();


            for(const auto& is : intersections(gridView_, directNeighbor))
            {
                constexpr int codimCommonEntity = 2;// TODO: 3d?

                auto isNormalSide = [&](const int idx)
                {
                    return !(idx == neighborIsIdx || localOppositeIdx_(idx) == neighborIsIdx);
                };

                auto localToGlobalEntityIdx = [&](const int localIdx)
                {
                    return gridView_.indexSet().subIndex(directNeighbor, localIdx, codimCommonEntity);
                };
                // skip the directly neighboring face itself
                if(isNormalSide(is.indexInInside()))
                {
                    // get facet
//                     const auto& facet = directNeighbor.template subEntity < 1 > (is.indexInInside());
//                     std::cout <<

                    // iterate over facets sub-entities // TODO: get number correctly
                    for(int i = 0; i < 2; ++i)
                    {
                        int localCommonEntIdx = referenceElement.subEntity(is.indexInInside(), 1, i, dim);
                        const auto& commonEnt = directNeighbor.template subEntity < codimCommonEntity > (localCommonEntIdx);

//                         int globalCommonEntIdx = gridView_.indexSet().subIndex(directNeighbor, localCommonEntIdx, codimCommonEntity);
                        int globalCommonEntIdx = localToGlobalEntityIdx(localCommonEntIdx);

//                         int dofIdx = gridView_.indexSet().index(commonEnt);

//                         if(globalCommonEntIdx != dofIdx)
//                             std::cout << "error!";

                        std::cout << "pos: " << commonEnt.geometry().center() << std::endl;

                        if(globalCommonEntIdx == localToGlobalEntityIdx(pairData_[0].localCommonEntIdx))
                        {
                            pairData_[0].normalPair.second = gridView_.indexSet().subIndex(directNeighbor, is.indexInInside(), dim-1) + offset_;
                        }
                        if(globalCommonEntIdx == localToGlobalEntityIdx(pairData_[1].localCommonEntIdx))
                        {
                            pairData_[1].normalPair.second = gridView_.indexSet().subIndex(directNeighbor, is.indexInInside(), dim-1) + offset_;
                        }
                    }

                    std::cout << "inters: " << intersection_.geometry().center() <<
                    " || side : " << is.geometry().center() << std::endl;


                }

            }
        }

//TODO: fill parallel faces
    }

private:

    int localOppositeIdx_(const int idx) const
    {
        return (idx % 2) ? (idx - 1) : (idx + 1);
    }


    template<class G = GridView>
    typename std::enable_if<G::dimension == 2, void>::type
    setInnerNormalPairs_(const int isIdx)
    {
        auto& data1 = pairData_[0];
        auto& data2 = pairData_[1];

        switch(isIdx)
        {
            case 0:
                data1.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 3, dim-1) + offset_;
                data2.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 2, dim-1) + offset_;
                data1.localCommonEntIdx = 2;
                data2.localCommonEntIdx = 0;
                break;
            case 1:
                data1.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 2, dim-1) + offset_;
                data2.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 3, dim-1) + offset_;
                data1.localCommonEntIdx = 1;
                data2.localCommonEntIdx = 3;
                break;
            case 2:
                data1.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 0, dim-1) + offset_;
                data2.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 1, dim-1) + offset_;
                data1.localCommonEntIdx = 0;
                data2.localCommonEntIdx = 1;
                break;
            case 3:
                data1.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 1, dim-1) + offset_;
                data2.normalPair.first = gridView_.indexSet().subIndex(intersection_.inside(), 0, dim-1) + offset_;
                data1.localCommonEntIdx = 3;
                data2.localCommonEntIdx = 2;
                break;
            default:
                DUNE_THROW(Dune::InvalidStateException, "foo");
        }
    }

    template<class G = GridView>
    typename std::enable_if<G::dimension == 3, void>::type
    setInnerNormalPairs_(const int isIdx)
    {
        // TODO: 3D
        DUNE_THROW(Dune::NotImplemented, "3d helper not ready yet");
    }


    const Intersection& intersection_;
    const typename Element::Geometry& elementGeometry_; //! Reference to the element geometry
    const GridView gridView_;
    const int offset_;

    std::array<PairData<Scalar>, numPairs> pairData_;


};

} // end namespace Dumux

#endif
