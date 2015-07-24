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
 * \brief Test for the coupled non-isothermal two-component ZeroEq and
 *        non-isothermal two-phase two-component Darcy model
 */

#include "config.h"
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_MULTIDOMAIN

#include <dumux/common/start.hh>

#include "2cnizeroeq2p2cniproblem.hh"

/*!
 * \brief Print a usage string for simulations.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void printUsage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
        errorMessageOut += progName;
        errorMessageOut += " [options]\n";
        errorMessageOut += errorMsg;
        errorMessageOut += "\nThe list of mandatory options for this program is:\n"
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
                           "NoDarcyX1                Horizontal position where the porous medium starts [m]\n"
                           "NoDarcyX2                Horizontal position where the porous medium ends [m]\n"
                           "RunUpDistanceX1          Horizontal position where the coupling starts [m]\n"
                           "RunUpDistanceX2          Horizontal position where the coupling ends [m]\n"
                           "\n"
                           "[SpatialParams]\n"
                           "AlphaBJ                  Beavers-Joseph coefficient [-]\n"
                           "Permeability             Hydraulic conductivity [m^2]\n"
                           "Porosity                 Porosity [-]\n"
                           "Swr                      Residual water saturation [-]\n"
                           "Snr                      Residual gas saturation [-]\n"
                           "VgAlpha                  Van-Genuchten parameter [1/Pa]\n"
                           "VgN                      Van-Genuchten parameter [-]\n"
                           "ThermalConductivitySolid Thermal conductivity of the solid material [W/(m*K)]\n"
                           "\n"
                           "[FreeFlow]\n"
                           "RefVelocity              Inflow velocity [m/s]\n"
                           "RefPressure              Reference pressure [Pa]\n"
                           "RefMassfrac              Inflow water mass fraction [-]\n"
                           "RefTemperature           Inflow temperature [K]\n"
                           "\n"
                           "[PorousMedium]\n"
                           "RefSw                    Initial water saturation [-]\n"
                           "RefPressurePM            Initial pressure [Pa]\n"
                           "RefTemperaturePM         Initial temperature [K]\n"
                           "\n"
                           "[Output]\n"
                           "NameFF                   Name free flow .vtu files\n"
                           "NamePM                   Name porous medium .vtu files\n"
                           "FreqRestart              Frequency of writting restart information\n"
                           "FreqOutput               Frequency of writting vtu output\n"
                           "FreqMassOutput           Frequency of writting storage output\n"
                           "FreqFluxOutput           Frequency of writting flux output\n"
                           "FreqVaporFluxOutput      Frequency of writting vapor flux output\n"
                           "\n"
                           "[TimeManager]\n"
                           "EpisodeLength            Length of one episode [s]\n"
                           "\n"
                           "[BoundaryLayer]\n"
                           "Model                    Enable use of boundary layer models (discouraged)\n"
                           "\n"
                           "[MassTransfer]\n"
                           "Model                    Enable use of mass transfer models (discouraged)\n"
                           "\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
#if (HAVE_SUPERLU || HAVE_UMFPACK)
    typedef TTAG(TwoCNIZeroEqTwoPTwoCNIProblem) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, printUsage);
#else
#warning "You need to have SuperLU or UMFPack installed to run this test."
    std::cerr << "You need to have SuperLU or UMFPack installed to run this test\n";
    return 77;
#endif
}

#else
int main(int argc, char** argv)
{
#warning You need to have dune-multidomain installed to run this test
    std::cerr << "You need to have dune-multidomain installed to run this test\n";
    return 77;
}
#endif
