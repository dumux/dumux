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
 *
 * \brief Creates a list of vertex indices on the process border which
 *        can be used to construct the foreign overlap.
 */
#ifndef DUMUX_VERTEX_BORDER_LIST_FROM_GRID_HH
#define DUMUX_VERTEX_BORDER_LIST_FROM_GRID_HH

#include "borderindex.hh"

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

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
template <class GridView, class VertexMapper>
class VertexBorderListFromGrid : public Dune::CommDataHandleIF<VertexBorderListFromGrid<GridView, VertexMapper>,
                                                               int >
{
    typedef std::list<BorderIndex> BorderList;

public:
    VertexBorderListFromGrid(const GridView &gridView,
                             const VertexMapper &map)
        : gridView_(gridView), map_(map)
    {
        gridView.communicate(*this,
                       Dune::InteriorBorder_InteriorBorder_Interface,
                       Dune::ForwardCommunication);
    };

    // data handle methods
    bool contains (int dim, int codim) const
    { return dim == codim; }

    bool fixedsize(int dim, int codim) const
    { return true; }

    template<class EntityType>
    size_t size(const EntityType &e) const
    { return 2; }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp &buff, const EntityType &e) const
    {
        buff.write(gridView_.comm().rank());
        buff.write(static_cast<int>(map_.map(e)));
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        BorderIndex bIdx;

        bIdx.localIdx = map_.map(e);
        buff.read(bIdx.peerRank);
        buff.read(bIdx.peerIdx);
        bIdx.borderDistance = 0;
        // vertices on the border are always in the interior of more
        // than one process which means that they are shared.
        bIdx.isShared = true;

        borderList_.push_back(bIdx);
    }

    // Access to the foreign border list.
    const BorderList &foreignBorderList() const
    { return borderList_; }

    // Access to the domestic border list (same as foreign border list
    // because all vertices are shared entities)
    const BorderList &domesticBorderList() const
    { return borderList_; }

private:
    const GridView gridView_;
    const VertexMapper &map_;
    BorderList borderList_;
};

} // namespace Dumux

#endif
