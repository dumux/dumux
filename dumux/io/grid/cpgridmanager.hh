// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief A grid creator that reads Petrel files and generates a CpGrid.
 */
#ifndef DUMUX_IO_GRID_CPGRIDMANAGER_HH
#define DUMUX_IO_GRID_CPGRIDMANAGER_HH

#include <config.h>

#if HAVE_OPM_GRID
#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/grid/CpGrid.hpp>

#if DUNE_VERSION_GTE(OPM_GRID, 2022, 10)
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#else
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <dumux/common/parameters.hh>

#if HAVE_ECL_INPUT

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
        Opm::EclipseState eclState(*deck_);
        grid_ = std::make_shared<Grid>();
        grid_->processEclipseFormat(&eclGrid, &eclState, false, false, false);
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
        if (Dune::MPIHelper::getCommunication().size() > 1)
            grid_->loadBalance();
    }

private:
    std::shared_ptr<Deck> deck_; //!< the eclipse deck
    std::shared_ptr<Grid> grid_; //!< the grid pointer
};

} // end namespace Dumux

#else
#warning "Eclipse input support in opm-common is required to use the cornerpoint grid manager"
#endif // HAVE_ECL_INPUT

#endif // HAVE_OPM_GRID

#endif
