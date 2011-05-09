// $Id$
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Creates a list of element indices on the process border which
 *        can be used to construct the foreign overlap.
 */
#ifndef DUMUX_ELEMENT_BORDER_LIST_FROM_GRID_HH
#define DUMUX_ELEMENT_BORDER_LIST_FROM_GRID_HH

#include <dune/grid/common/datahandleif.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include "borderindex.hh"

#include <algorithm>
#include <list>
#include <set>
#include <map>

namespace Dumux {


/*!
 * \brief Uses communication on the grid to find the initial seed list
 *        of indices.
 *
 * \todo implement this class generically. For this, it must be
 *       possible to query the mapper whether it contains entities of
 *       a given codimension without the need to hand it an actual
 *       entity.
 */
template <class GridView, class ElementMapper>
class ElementBorderListFromGrid : public Dune::CommDataHandleIF<ElementBorderListFromGrid<GridView, ElementMapper>,
                                                                int >
{
    typedef int ProcessRank;
    typedef int Index;
    typedef int Distance;
    typedef Index LocalIndex;
    typedef Index PeerIndex;
    typedef std::list<BorderIndex> BorderList;

public:
    ElementBorderListFromGrid(const GridView &gv,
                             const ElementMapper &map)
        : gv_(gv), map_(map)
    {
        gv.communicate(*this, 
                       Dune::InteriorBorder_All_Interface,
//                       Dune::Overlap_All_Interface,
                       Dune::BackwardCommunication);
    };

    // data handle methods
    bool contains (int dim, int codim) const
    { return codim == 0; }
    
    bool fixedsize(int dim, int codim) const
    { return true; }
    
    template<class EntityType> 
    size_t size(const EntityType &e) const
    { return 2; };

    template<class MessageBufferImp, class EntityType> 
    void gather(MessageBufferImp &buff, const EntityType &e) const 
    { 
        buff.write(gv_.comm().rank()); 
        buff.write(map_.map(e)); 
    };

    template<class MessageBufferImp, class EntityType> 
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        BorderIndex bIdx;
        
        bIdx.localIdx = map_.map(e);
        buff.read(bIdx.peerRank);
        buff.read(bIdx.peerIdx);
        bIdx.borderDistance = 1;

        borderList_.push_back(bIdx);
    };
    
    // Access to the initial seed list.
    const BorderList &borderList() const
    { return borderList_; }

private:
    const GridView gv_;
    const ElementMapper &map_;
    BorderList borderList_;
};

} // namespace Dumux

#endif
