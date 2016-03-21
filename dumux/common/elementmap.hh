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
 * \brief Base class for the finite volume geometry vector for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
#ifndef DUMUX_ELEMENT_MAP_HH
#define DUMUX_ELEMENT_MAP_HH

namespace Dumux {

//! An index to element map
template <class GridView>
class ElementMap
  : public std::vector<typename GridView::Traits::Grid::template Codim<0>::EntitySeed>
{
    using Grid = typename GridView::Traits::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
public:
    ElementMap(const GridView& gridView_) : grid_(gridView_.grid()) {}

    Element element(IndexType eIdx) const
    { return grid_.entity((*this)[eIdx]); }

private:
    const Grid& grid_;
};

} // end namespace Dumux

#endif