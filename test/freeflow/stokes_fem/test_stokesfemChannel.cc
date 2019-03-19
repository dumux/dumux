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
 * \brief Test for the isothermal Stokes box model.
 */
#include <config.h>
#include "stokesfemtestproblemChannel.hh"
#include <dumux/common/start.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe List of Mandatory arguments for this program is:\n"
                                        "\t-TimeManager.TEnd               The end of the simulation. [s] \n"
                                        "\t-TimeManager.DtInitial          The initial timestep size. [s] \n"
                                        "\t-Grid.File                      The file name of the file containing the grid \n"
                                        "\t                                   definition in DGF format\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
#if (HAVE_SUPERLU || HAVE_UMFPACK)
    typedef TTAG(StokesFemTestChannel) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv, usage);
#else
#warning "You need to have SuperLU or UMFPack installed to run this test."
    std::cerr << "You need to have SuperLU or UMFPack installed to run this test\n";
    return 77;
#endif
}
