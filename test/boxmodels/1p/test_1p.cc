// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 *
 * \brief test for the one-phase box model
 */
#include "config.h"
#include "1ptestproblem.hh"
#include <dumux/common/start.hh>

void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::cout << errorMsg << "\n"
                  << "\n";
    }
    std::cout
        << "Usage: " << progName << " [options]\n"
        << "Mandatory options are:\n"
        << "\t--t-end=ENDTIME                  The time of the end of the simlation [s]\n"
        << "\t--dt-initial=STEPSIZE            The initial time step size [s]\n"
        << "\t--grid-file=FILENAME             The file name of the file containing the grid \n"
        << "\t                                   definition in DGF format\n"
        << "\n"
        << "Alternativ supported syntax:\n"
        << "\t-tEnd ENDTIME                    The time of the end of the simlation [s]\n"
        << "\t-dtInitial STEPSIZE              The initial time step size [s]\n"
        << "\t-gridFile FILENAME              The file name of the file containing the grid \n"
        << "\t                                   definition in DGF format\n"
        << "\n"
        << "If --parameter-file is specified parameters can also be defined there. In this case,\n"
        << "camel case is used for the parameters (e.g.: --grid-file becomes gridFile). Parameters\n"
        << "specified on the command line have priority over those in the parameter file.\n"
        << "Important optional options include:\n"
        << "\t--help,-h                        Print this usage message and exit\n"
        << "\t--print-parameters[=true|false]  Print the run-time modifiable parameters _after_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--print-properties[=true|false]  Print the compile-time parameters _before_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--parameter-file=FILENAME        File with parameter definitions\n"
        << "\t--restart=RESTARTTIME            Restart simulation from a restart file\n"
        << "\n"
        << "For the case of no arguments given, the input parameter file is expected to be named './parameter.input' \n"
        << "\n";
}

int main(int argc, char** argv)
{
    typedef TTAG(OnePTestProblem) ProblemTypeTag;
    return Dumux::startWithParameters<ProblemTypeTag>(argc, argv, usage);
}
