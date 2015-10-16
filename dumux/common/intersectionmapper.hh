// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
#ifndef DUMUX_INTERSECTIONITERATOR_HH
#define DUMUX_INTERSECTIONITERATOR_HH

#include<vector>
#include<unordered_map>
#include<dune/grid/common/mcmgmapper.hh>

/**
 * @file
 * @brief  defines an intersection mapper for mapping of global DOF's assigned to faces which also works for adaptive grids.
 */

namespace Dumux
{

template<class GridView>
class IntersectionMapper
{
    typedef typename GridView::Grid Grid;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;

public:
    IntersectionMapper(const GridView& gridview)
    : gridView_(gridview), elementMapper_(gridView_), size_(gridView_.size(1)),
      intersectionMapGlobal_(gridView_.size(0)), intersectionMapLocal_(gridView_.size(0))
    {
        const auto element = *gridView_.template begin<0>();

        int fIdx = 0;
        for (const auto& intersection : Dune::intersections(gridView_, element))
        {
            int idxInInside = intersection.indexInInside();

            standardLocalIdxMap_[idxInInside] = fIdx;

            fIdx++;
        }
    }

    const ElementMapper& elementMapper() const
    {
        return elementMapper_;
    }

    int map(const Element& element) const
    {
        return elementMapper_.index(element);
    }

    int map(int elemIdx, int fIdx)
    {
        return intersectionMapGlobal_[elemIdx][fIdx];
    }

    int map(int elemIdx, int fIdx) const
    {
        return (intersectionMapGlobal_[elemIdx].find(fIdx))->second;//use find() for const function!
    }

    int map(const Element& element, int fIdx)
    {
        return intersectionMapGlobal_[map(element)][fIdx];
    }

    int map(const Element& element, int fIdx) const
    {
        return intersectionMapGlobal_[map(element)].find(fIdx)->second;//use find() for const function!
    }

    int maplocal(int elemIdx, int fIdx)
    {
        return intersectionMapLocal_[elemIdx][fIdx];
    }

    int maplocal(int elemIdx, int fIdx) const
    {
        return (intersectionMapLocal_[elemIdx].find(fIdx))->second;//use find() for const function!
    }


    int maplocal(const Element& element, int fIdx)
    {
        return intersectionMapLocal_[map(element)][fIdx];
    }

    int maplocal(const Element& element, int fIdx) const
    {
        return (intersectionMapLocal_[map(element)].find(fIdx))->second;//use find() for const function!
    }

    // return number intersections
    unsigned int size () const
    {
        return size_;
    }

    unsigned int size (int elemIdx) const
    {
        return intersectionMapLocal_[elemIdx].size();
    }

    unsigned int size (const Element& element) const
    {
        return intersectionMapLocal_[map(element)].size();
    }

    void update()
    {
        elementMapper_.update();

        intersectionMapGlobal_.clear();
        intersectionMapGlobal_.resize(elementMapper_.size());
        intersectionMapLocal_.clear();
        intersectionMapLocal_.resize(elementMapper_.size());

        for (const auto& element : Dune::elements(gridView_))
        {
            int eIdxGlobal = map(element);

            int fIdx = 0;
            // run through all intersections with neighbors
            for (const auto& intersection : Dune::intersections(gridView_, element))
            {
                int indexInInside = intersection.indexInInside();
                intersectionMapLocal_[eIdxGlobal][fIdx] = indexInInside;

                fIdx++;
            }
        }

        int globalIntersectionIdx = 0;
        for (const auto& element : Dune::elements(gridView_))
        {
            int eIdxGlobal = map(element);

            int fIdx = 0;
            // run through all intersections with neighbors
            for (const auto& intersection : Dune::intersections(gridView_, element))
            {
                if (intersection.neighbor())
                {
                    auto neighbor = intersection.outside();
                    int globalIdxNeighbor = map(neighbor);

                    if (element.level() > neighbor.level() || (element.level() == neighbor.level() && eIdxGlobal < globalIdxNeighbor))
                    {

                        int faceIdxNeighbor = 0;
                        if (size(globalIdxNeighbor) == 2 * dim)
                        {
                            faceIdxNeighbor = standardLocalIdxMap_[intersection.indexInOutside()];
                        }
                        else
                        {
                            for (const auto& intersectionNeighbor : Dune::intersections(gridView_, neighbor))
                            {
                                if (intersectionNeighbor.neighbor())
                                {
                                    if (intersectionNeighbor.outside() == element)
                                    {
                                        break;
                                    }
                                }
                                faceIdxNeighbor++;
                            }
                        }

                        intersectionMapGlobal_[eIdxGlobal][fIdx] = globalIntersectionIdx;
                        intersectionMapGlobal_[globalIdxNeighbor][faceIdxNeighbor] = globalIntersectionIdx;
                        globalIntersectionIdx ++;
                    }
                }
                else
                {
                    intersectionMapGlobal_[eIdxGlobal][fIdx] = globalIntersectionIdx;
                    globalIntersectionIdx ++;
                }
                fIdx++;
            }
        }
        size_ = globalIntersectionIdx;
    }

protected:
    const GridView& gridView_;
    ElementMapper elementMapper_;
    unsigned int size_;
    std::vector<std::unordered_map<int, int> > intersectionMapGlobal_;
    std::vector<std::unordered_map<int, int> > intersectionMapLocal_;
    std::unordered_map<int, int> standardLocalIdxMap_;
};

}

#endif
