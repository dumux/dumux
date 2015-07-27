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
/**
 * \file
 *
 * \brief Test for the kinetic modules of the mpnc box model.
 */
#include "config.h"

#include <dumux/common/start.hh>
#include "evaporationatmosphereproblem.hh"


/*!
 * \brief Print a usage string for simulations.
 *
 * \param progName The name of the program, that was tried to be started.
 * \param errorMsg The error message that was issued by the start function.
 *                 Comprises the thing that went wrong and a general help message.
 */
void printUsage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
        errorMessageOut += progName;
        errorMessageOut += " [options]\n";
        errorMessageOut += errorMsg;
        errorMessageOut += "\nAn uncomplete list of mandatory options for this program is:\n"
                           "[Grid]\n"
                           "LowerLeftX               Minumum x-coordinate [m]\n"
                           "UpperRightX              Maximum x-coordinate [m]\n"
                           "LowerLeftY               Minumum y-coordinate [m]\n"
                           "UpperRightY              Maximum y-coordinate [m]\n"
                           "NumberOfCellsX           Number of cells in x-direction\n"
                           "NumberOfCellsY           Number of cells in y-direction\n"
                           "GradingFactorY           Vertical grading of the cells\n"
                           "RefineTop                Specifies whethter the top of the free flow will be refined\n"
                           "InterfacePosY            Vertical position of the interface [m]\n"
                           "\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
#if !HAVE_UG && !(HAVE_ALUGRID || HAVE_DUNE_ALUGRID)
#warning Evaporation Atmosphere not built, needs either UG or dune-ALUGrid for the log mesh.
    std::cerr << "Evaporation Atmosphere not built, needs either UG or dune-ALUGrid for the log mesh." << std::endl;
    return 77;
#else
    typedef TTAG(EvaporationAtmosphereProblem) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, printUsage);//, usage);
#endif
}
