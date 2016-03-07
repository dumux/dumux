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
 * \brief Main file of the tutorial for a fully coupled twophase box model.
 */
#include <config.h> /*@\label{tutorial-implicit:include-begin}@*/
#include "tutorialproblem_implicit.hh"  /*@\label{tutorial-implicit:include-problem-header}@*/
#include <dumux/common/start.hh> /*@\label{tutorial-implicit:include-end}@*/

//! Prints a usage/help message if something goes wrong or the user asks for help
void usage(const char *progName, const std::string &errorMsg)  /*@\label{tutorial-implicit:usage-function}@*/
{
    std::cout
        <<  "\nUsage: " << progName << " [options]\n";
    if (errorMsg.size() > 0)
        std::cout << errorMsg << "\n";
    std::cout
        << "\n"
        << "The list of mandatory arguments for this program is:\n"
        << "\t-TEnd                The end of the simulation [s]\n"
        << "\t-DtInitial           The initial timestep size [s]\n"
        << "\t-Grid.UpperRightX    The x-coordinate of the grid's upper-right corner [m]\n"
        << "\t-Grid.UpperRightY    The y-coordinate of the grid's upper-right corner [m]\n"
        << "\t-Grid.NumberOfCellsX The grid's x-resolution\n"
        << "\t-Grid.NumberOfCellsY The grid's y-resolution\n"
        << "\n";
}

int main(int argc, char** argv)
{
    typedef TTAG(TutorialProblemImplicit) TypeTag; /*@\label{tutorial-implicit:set-type-tag}@*/
    return Dumux::start<TypeTag>(argc, argv, usage); /*@\label{tutorial-implicit:call-start}@*/
}
