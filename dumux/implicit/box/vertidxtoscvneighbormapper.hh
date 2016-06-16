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
#ifndef DUMUX_VertIdxToScvNeighborMapper_HH
#define DUMUX_VertIdxToScvNeighborMapper_HH

#include<vector>
#include<unordered_map>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/grid/common/grid.hh>
#include <fstream>
using namespace std;

/**
 * @file
 * @brief  defines a vertex mapper for mapping neighbor elements and neighbor subcontrol volumes to a global vertex index.
 */

namespace Dumux
{

template<class GridView>
class VertIdxToScvNeighborMapper
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
    // iterator for the neighbor element pointer and scvIdx pair
    typedef typename NeighborElementScv::const_iterator NeighborElementScvIterator;
    // map containing vertex index as a key and a vector of pairs of element pointers and neighbor scv index
    typedef std::unordered_map<int, NeighborElementScv > VertexElementScvMapper;

    VertIdxToScvNeighborMapper(const GridView& gridView)
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
    unsigned int vertexElementsScvIdx (int vIdxGlobal, int elem) const
    {
        return vertexElements(vIdxGlobal)[elem].second;
    }

    //return pointer to i-th element
    ElementPointer vertexElementPointer (int vIdxGlobal, int elem) const
    {
        return gridView_.grid().entity(vertexElements(vIdxGlobal)[elem].first);
    }

     /*!
     * \brief Updates the VertexToScvNeighborMapper
     * \note  Loops over Elements to match subcontrol volumes to the respective adjacent vertices
     */
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
                int vIdxGlobal = vertexMapper_.subIndex(*eIt, vIdx, dim);
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

    /*!
    * \brief Produces a text file which lists for each vertex the neighbored elements as well as the indices of
    * the vertex adjacent subcontrol volumes
    */
    void output()
    {
        vertexMapper_.update();
        int numDofs = gridView_.size(dim);
        ofstream out("vertexScvMapper.txt", ios::app);
        for (int globalIdx = 0; globalIdx < numDofs; ++globalIdx)
            {
                out << "global vertex index: " << globalIdx << endl;
                int numNeighbors = size(globalIdx);
                out << "number of neighbor elements: " << numNeighbors << endl;
                for (int neighborIdx = 0; neighborIdx < numNeighbors; neighborIdx++) {
                    int neighborScvIdx = vertexElementsScvIdx(globalIdx, neighborIdx);
                    ElementPointer elem = vertexElementPointer(globalIdx, neighborIdx);
                    out << neighborIdx << "-> element center: (" << elem.geometry().center()[0] << ", " << elem.geometry().center()[1] << ") ";
                    out << "scvIdx: " << neighborScvIdx << endl;
                }
                out << "---------------" << endl;
            }
        out.close();
    }

protected:

    const GridView& gridView_;
    VertexElementScvMapper vertexElementScvMapper_;
    //vertex map
    VertexMapper vertexMapper_;

};

}

#endif
