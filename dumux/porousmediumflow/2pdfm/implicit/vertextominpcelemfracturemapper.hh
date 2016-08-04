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
#ifndef DUMUX_VertIdxToMinPcFractureMapper_HH
#define DUMUX_VertIdxToMinPcFractureMapper_HH

#include<vector>
#include<unordered_map>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/grid/common/grid.hh>
#include <fstream>
#include "vertidxtoelemfacesmapper.hh"
using namespace std;

/**
 * @file
 * @brief  defines a vertex mapper for mapping neighbor elements and neighbor subcontrol volumes to a global vertex index.
 */

namespace Dumux
{

template<class TypeTag>
class VertIdxToMinPcFractureMapper
{

typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::template Codim<0>::Entity::EntitySeed EntitySeed;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGVertexLayout> VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Dumux::VertIdxToElemFacesMapper<GridView> VertIdxToElemFacesMapper;



public:

    //public types
    // vector of pairs of element seed and minPe
    typedef std::pair<EntitySeed, Scalar>  ElementMinPe;
    // map containing vertex index as a key and a vector of pairs of element pointers and neighbor scv index
    typedef std::unordered_map<int, ElementMinPe > VertexElemMinPcFractureMapper;

    VertIdxToMinPcFractureMapper(const GridView& gridView, SpatialParams& spatialParams)
    : gridView_(gridView), vertexMapper_(gridView), spatialParams_(spatialParams), vertIdxToElemFacesMapper_(gridView)

    {
        //Call the update function in the constructor
    }

    // reference to the neighbor elements scvIcx pairs
    const ElementMinPe& elementMinPe(int vIdxGlobal) const
    {
        return vertexElemMinPcFractureMapper_.find(vIdxGlobal)->second;
    }

    // return the reference to the mapper
    const VertexElemMinPcFractureMapper& map() const
    {
        return vertexElemMinPcFractureMapper_;
    }

    // return map size (number of nodes)
    unsigned int size () const
    {
        return vertexElemMinPcFractureMapper_.size();
    }

    //return current Pe for given Vertex
    Scalar currentMinPe (int vIdxGlobal) const
    {
        return elementMinPe(vIdxGlobal).second;
    }

    //return pointer to i-th element
    Element vertexElementPointer (int vIdxGlobal) const
    {
        return gridView_.grid().entity(elementMinPe(vIdxGlobal).first);
    }

     /*!
     * \brief Updates the VertexToScvNeighborMapper
     * \note  Loops over Elements to match subcontrol volumes to the respective adjacent vertices
     */
    void update()
    {

        updateBrooksCorey();
    }


    void updateBrooksCorey()
    {
        vertexMapper_.update();
        vertIdxToElemFacesMapper_.update();
        vertexElemMinPcFractureMapper_.clear();
        FVElementGeometry fvGeometry;
        ElementIterator eEndIt = gridView_.template end<0>();
        //Loop over elements
        for (ElementIterator eIt = gridView_.template begin<0>(); eIt != eEndIt; ++eIt)
        {
            //get the number of vertices in the element
            int noCorner = eIt->geometry().corners();
            //create and store the entity seed which uses as little memory as possible
            EntitySeed entitySeed = eIt->seed();
            Scalar pe;
            //Loop over local nodes (respectively scvIdx) in element
            fvGeometry.update(gridView_, *eIt);
            for (int vIdx = 0; vIdx < noCorner; ++vIdx)
            {
                //Get the global vertex index
                int vIdxGlobal = vertexMapper_.subIndex(*eIt, vIdx, dim);
                typename VertexElemMinPcFractureMapper::iterator vHandle = vertexElemMinPcFractureMapper_.find(vIdxGlobal);

                bool isFracture = spatialParams_.isVertexFracture(*eIt, vIdx);
                if (isFracture)
                    pe = spatialParams_.materialLawParamsFracture(*eIt, fvGeometry, vIdx).pe();
                else
                    pe = spatialParams_.materialLawParams(*eIt, fvGeometry, vIdx).pe();
                //   if vertex (number) was not found yet, insert new vertex entry with first element neighbor
                if (vHandle == vertexElemMinPcFractureMapper_.end())
                {
                    vertexElemMinPcFractureMapper_.insert(std::make_pair(vIdxGlobal, (std::make_pair(entitySeed, pe))));
                }
                // else if vertex index was already found, insert another element neighbor to that vertex
                else if (pe < currentMinPe(vIdxGlobal))
                {
                    vertexElemMinPcFractureMapper_[vIdxGlobal] = std::make_pair(entitySeed, pe);
                }
            }
        }
    }

  //  void updatevanGenuchten()
  //  {
  //      vertexMapper_.update();
  //      vertIdxToElemFacesMapper_.update();
  //      vertexElemMinPcFractureMapper_.clear();
  //      FVElementGeometry fvGeometry;
  //      ElementIterator eEndIt = gridView_.template end<0>();
  //      //Loop over elements
  //      for (ElementIterator eIt = gridView_.template begin<0>(); eIt != eEndIt; ++eIt)
  //      {
  //          int noCorner = eIt->geometry().corners();
  //          //create and store the entity seed which uses as little memory as possible
  //          EntitySeed entitySeed = eIt->seed();
  //          //Loop over local nodes (respectively scvIdx) in element
  //          fvGeometry.update(gridView_, *eIt);
  //          //for (int vIdx = 0; vIdx < noCorner; ++vIdx)
  //          for (int vIdx = noCorner-1; vIdx >= 0; vIdx--)
  //          {
  //              cout << vIdx << std::endl;
  //              //Get the global vertex index
  //              int vIdxGlobal = vertexMapper_.subIndex(*eIt, vIdx, dim);
  //              // check if there is an adjacent fracture
  //              bool isFracture = spatialParams_.isVertexFracture(*eIt, vIdx);
  //              typename VertexElemMinPcFractureMapper::iterator vHandle = vertexElemMinPcFractureMapper_.find(vIdxGlobal);
  //              //   if vertex (number) was not found yet, insert new vertex entry with first element neighbor
  //              if (vHandle == vertexElemMinPcFractureMapper_.end())
  //              {
  //                  vertexElemMinPcFractureMapper_.insert(std::make_pair(vIdxGlobal, (std::make_pair(entitySeed, 0))));
  //              }
  //          }
  //      }
  //  }


    /*!
    * \brief Produces a text file which lists for each vertex the neighbored elements as well as the indices of
    * the vertex adjacent subcontrol volumes
    */
    void output()
    {
        vertexMapper_.update();
        int numDofs = gridView_.size(dim);
        ofstream out("vertexMinPc.txt", ios::app);
        for (int globalIdx = 0; globalIdx < numDofs; ++globalIdx)
            {
                out << "global vertex index: " << globalIdx << endl;
                    out << "currentMinPe" << currentMinPe(globalIdx) << endl;
                out << "---------------" << endl;
            }
        out.close();
    }

protected:

    const GridView& gridView_;
    VertexElemMinPcFractureMapper vertexElemMinPcFractureMapper_;
    //vertex map
    VertexMapper vertexMapper_;
    SpatialParams spatialParams_;
    VertIdxToElemFacesMapper vertIdxToElemFacesMapper_;

};

}

#endif
