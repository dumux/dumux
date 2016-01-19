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
 * \brief Test for the coupled isothermal two-component Stokes and
 *        isothermal two-phase two-component Darcy model
 */

#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_MULTIDOMAIN

#include <dumux/common/start.hh>

#include "2cstokes2p2cproblem.hh"

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
        errorMessageOut += "\n\nThe list of optional options for this program is:\n"
                           "[BoundaryLayer]\n"
                           "Model               Number/ID of the used model\n"
                           "Offset              Virtual run-up distance for BL models [m]\n"
                           "ConstThickness      Constant BL thickness (model 1) [m]\n"
                           "YPlus               Conversion value (model 4-6) [-]\n"
                           "RoughnessLength     Characteristic roughness length (model 6)\n"
                           "\n"
                           "[MassTransferModel]\n"
                           "Coefficient         Coeffient used for the exponential law (model 1) [-]\n"
                           "CharPoreRadius      Characteristic pore radius for Schluender model (model 2+4) [m]\n"
                           "\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
#if (HAVE_SUPERLU || HAVE_PARDISO)
    typedef TTAG(TwoCStokesTwoPTwoCProblem) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, printUsage);
#else
#warning "You need to have SuperLU or Pardiso installed to run this test."
    std::cerr << "You need to have SuperLU or Pardiso installed to run this test\n";
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
