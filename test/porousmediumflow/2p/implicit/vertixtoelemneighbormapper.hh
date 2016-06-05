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
#ifndef DUMUX_VERTIDXTOELEMNEIGHBORMAPPER_HH
#define DUMUX_VERTIDXTOELEMNEIGHBORMAPPER_HH

#include<vector>
#include<unordered_map>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/grid/common/grid.hh>

/**
 * @file
 * @brief  defines a vertex mapper for mapping neighbor elements to a global vertex index. Every element in which a vertex appears is a neighbor element.
 */

namespace Dumux
{

template<class GridView>
class VertIdxToElemNeighborMapper
{

    typedef typename GridView::Grid Grid;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> VertexMapper;
    typedef typename Grid::template Codim<0>::Entity::EntitySeed EntitySeed;


public:

    //public types
    // vector of pairs of element seed and scvIdx
    typedef std::vector<std::pair<EntitySeed, int> > NeighborElementScv;
    // vector of element seeds
    //typedef std::vector<EntitySeed> NeighborElementSeeds;
    // iterator for the neighbor element pointer
    //typedef typename NeighborElementSeeds::const_iterator NeighborElementSeedsIterator;
    // iterator for the neighbor element pointer and scvIdx pair
    typedef typename NeighborElementScv::const_iterator NeighborElementScvIterator;
    // map containing vertex index as key and a vector neighbor element pointers as values
    //typedef std::unordered_map<int, NeighborElementSeeds > VertexElementMapper;
    // map for vertex index to pair of element pointers and scv index
    typedef std::unordered_map<int, NeighborElementScv > VertexElementScvMapper;

    VertIdxToElemNeighborMapper(const GridView& gridView)
    : gridView_(gridView), vertexMapper_(gridView)

    {
        //Call the update function in the constructor
        update();
    }

    // reference to the neighbor elements scvIcx pairs
    const NeighborElementScv& vertexElements(int vIdxGlobal) const
    {
        return vertexElementScvMapper_.find(vIdxGlobal)->second;
    }

    //return begin iterator
    NeighborElementScvIterator vertexElementsBegin(int vIdxGlobal) const
    {
        return vertexElementScvMapper_.find(vIdxGlobal)->second.begin();
    }

    //return end iterator
    NeighborElementScvIterator vertexElementsEnd(int vIdxGlobal) const
    {
        return vertexElementScvMapper_.find(vIdxGlobal)->second.end();
    }

    // return the reference to the mapper
    const VertexElementScvMapper& map() const
    {
        return vertexElementScvMapper_;
    }

    // return map size (number of nodes)
    unsigned int size () const
    {
        return vertexElementScvMapper_.size();
    }

    // return number of neighbor elements for a global vertex
    unsigned int size (int vIdxGlobal) const
    {
        return vertexElementScvMapper_.find(vIdxGlobal)->second.size();
    }

    //return scvIdx corresponding to i-th element
    unsigned int vertexElementsScvIdx (int vIdxGlobal, int elem)
    {
        return vertexElements(vIdxGlobal)[elem].second;
    }

    //return pointer to i-th element
    ElementPointer vertexElementPointer (int vIdxGlobal, int elem)
    {
        return gridView_.grid().entity(vertexElements(vIdxGlobal)[elem].first);
    }


    void update()
    {
        vertexMapper_.update();
        vertexElementScvMapper_.clear();

        ElementIterator eEndIt = gridView_.template end<0>();
        //Loop over elements
        for (ElementIterator eIt = gridView_.template begin<0>(); eIt != eEndIt; ++eIt)
        {
            //get the number of vertices in the element
            int noCorner = eIt->geometry().corners();
            //create and store the entity seed which uses as little memory as possible
            EntitySeed entitySeed = eIt->seed();
            //Loop over local nodes (respectively scvIdx) in element
            for (int vIdx = 0; vIdx < noCorner; ++vIdx)
            {
                //Get the global vertex index
//#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                int vIdxGlobal = vertexMapper_.subIndex(*eIt, vIdx, dim);

//#else
//                int vIdxGlobal = vertexMapper_.map(*eIt, vIdx, dim);
//#endif
                typename VertexElementScvMapper::iterator vHandle = vertexElementScvMapper_.find(vIdxGlobal);
                //   if vertex (number) was not found yet, insert new vertex entry with first element neighbor
                if (vHandle == vertexElementScvMapper_.end())
                {
                    vertexElementScvMapper_.insert(std::make_pair(vIdxGlobal, NeighborElementScv(1, std::make_pair(entitySeed, vIdx))));
                }
                // else if vertex index was already found, insert another element neighbor to that vertex
                else
                {
                    vHandle->second.push_back(std::make_pair(entitySeed, vIdx));
                }
            }
        }
    }

protected:

    const GridView& gridView_;
    //VertexElementMapper vertexElementMapper_;
    VertexElementScvMapper vertexElementScvMapper_;
    //vertex map
    VertexMapper vertexMapper_;

};

}

#endif
