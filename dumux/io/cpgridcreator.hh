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
 * \ingroup InputOutput
 * \brief A grid creator that reads Petrel files and generates a CpGrid.
 */
#ifndef DUMUX_CPGRID_CREATOR_HH
#define DUMUX_CPGRID_CREATOR_HH

#if HAVE_OPM_GRID
#include <dune/grid/CpGrid.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
}

/*!
 * \ingroup InputOutput
 * \brief A grid creator that reads Petrel files and generates a CpGrid.
 */
template <class TypeTag>
class CpGridCreator
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridPointer = std::shared_ptr<Grid>;
    using Deck = Opm::Deck;

public:
    /*!
     * \brief Create the Grid.
     */
    static void makeGrid()
    {
        std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Grid, File);

        deck() = Opm::Parser().parseFile(fileName);
        Opm::EclipseGrid ecl_grid(deck());

        gridPtr() = std::make_shared<Grid>(*(new Grid()));
        gridPtr()->processEclipseFormat(ecl_grid, false, false);
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
        static GridPointer cpGrid;
        return cpGrid;
    }

    /*!
     * \brief Returns a reference to the input deck.
     *
     * The input deck can be used to read parameters like porosity/permeability.
     */
    static Deck &deck()
    {
        static Deck deck_;
        return deck_;
    }

    /*!
     * \brief Distributes the grid over all processes for a parallel computation.
     */
    static void loadBalance()
    {
        gridPtr()->loadBalance();
    }
};
}
#endif // HAVE_OPM_GRID

#endif
