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
 * \brief Provides a grid creator which reads Dune Grid Format (DGF) files
 */
#ifndef DUMUX_DGF_GRID_CREATOR_HH
#define DUMUX_DGF_GRID_CREATOR_HH

#warning This file is deprecated and will be removed after Dumux 2.9. Use the default grid creator (Dumux::GridCreator).

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(Grid);
}

/*!
 * \brief Provides a grid creator which reads Dune Grid Format (DGF) files
 */
template <class TypeTag>
class DgfGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef Dune::GridPtr<Grid> GridPointer;

public:
    /*!
     * \brief Load the grid from the dgf file.
     */
    static void makeGrid(const std::string& dgfFileName)
    {
        gridPtr() = GridPointer(dgfFileName.c_str(), Dune::MPIHelper::getCommunicator());
    }

    /*!
     * \brief Load the grid from the dgf file given in the input file.
     */
    static void makeGrid()
    {
        const std::string dgfFileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Grid, File);
        makeGrid(dgfFileName);
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *gridPtr();
    }

    /*!
     * \brief Returns a reference to the grid pointer.
     */
    static GridPointer &gridPtr()
    {
        static GridPointer gridPtr_;
        return gridPtr_;
    }

    /*!
     * \brief Call loadBalance() function of GridPointer.
     */
    static void loadBalance()
    {
        gridPtr().loadBalance();
    }
};

} // namespace Dumux

#endif
