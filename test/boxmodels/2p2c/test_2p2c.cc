// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief test for the 2p2c box model
 */
#include "config.h"
#include "injectionproblem.hh"
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
        << "\t--dgf-file=FILENAME              The file name of the file containing the grid \n"
        << "\t                                 definition in DGF format\n"
        << "\n"
        << "Important optional options include:\n"
        << "\t--help,-h                        Print this usage message and exit\n"
        << "\t--print-parameters[=true|false]  Print the run-time modifiable parameters _after_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--print-properties[=true|false]  Print the compile-time parameters _before_ \n"
        << "\t                                 the simulation [default: true]\n"
        << "\t--opts-file=FILENAME             File with parameter definitions\n"
        << "\t--restart=RESTARTTIME            Restart simulation from a restart file\n"
        << "\n"
        << "If --opts-file is specified parameters can also be defined there. In this case,\n"
        << "camel case is used for the parameters (e.g.: --dgf-file becomes DgfFile). Parameters\n"
        << "specified on the command line have priority over those in the option file.\n" 
        << "\n";
}

int main(int argc, char** argv)
{
    typedef TTAG(InjectionProblem) ProblemTypeTag;
    return Dumux::startWithParameters<ProblemTypeTag>(argc, argv, usage);
}
