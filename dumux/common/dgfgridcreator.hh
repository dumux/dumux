// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Provides a grid creator which reads Dune Grid Format (DGF) files
 */
#ifndef DUMUX_DGF_GRID_CREATOR_HH
#define DUMUX_DGF_GRID_CREATOR_HH

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
     * \brief Load the grid from the file.
     */
    static void makeGrid()
    {
        const std::string dgfFileName = GET_RUNTIME_PARAM(TypeTag, std::string, DgfFile);

        gridPtr_ = GridPointer(dgfFileName.c_str());
    };

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *gridPtr_;
    };

    /*!
     * \brief Returns a reference to the grid pointer.
     *
     * This method is specific to the DgfGridCreator!
     */
    static GridPointer &gridPtr()
    {
        return gridPtr_;
    };
    
private:
    static GridPointer gridPtr_;
};

template <class TypeTag>
typename DgfGridCreator<TypeTag>::GridPointer DgfGridCreator<TypeTag>::gridPtr_;

} // namespace Dumux

#endif
