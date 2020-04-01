// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
#ifndef DUMUX_IO_GRID_CPGRIDMANAGER_HH
#define DUMUX_IO_GRID_CPGRIDMANAGER_HH

#if HAVE_OPM_GRID
#include <dune/common/parallel/mpihelper.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A grid creator that reads Petrel files and generates a CpGrid.
 */
class CpGridManager
{
public:
    using Grid = Dune::CpGrid;
    using Deck = Opm::Deck;

    /*!
     * \brief Create the Grid.
     */
    void init(const std::string& paramGroup = "")
    {
        const auto fileName = getParamFromGroup<std::string>(paramGroup, "Grid.File");
        deck_ = std::make_shared<Opm::Deck>(Opm::Parser().parseFile(fileName));
        Opm::EclipseGrid eclGrid(*deck_);
        grid_ = std::make_shared<Grid>();
        grid_->processEclipseFormat(&eclGrid, false, false, false);
        loadBalance();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid &grid()
    {
        return *grid_;
    }

    /*!
     * \brief Returns a reference to the input deck.
     *
     * The input deck can be used to read parameters like porosity/permeability.
     */
    std::shared_ptr<Deck> getDeck() const
    {
        return deck_;
    }

    /*!
     * \brief Distributes the grid over all processes for a parallel computation.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
            grid_->loadBalance();
    }

private:
    std::shared_ptr<Deck> deck_; //!< the eclipse deck
    std::shared_ptr<Grid> grid_; //!< the grid pointer
};

} // end namespace Dumux

#endif // HAVE_OPM_GRID

#endif
