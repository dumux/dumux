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
 * \brief test for the two-phase box model
 */
#include <config.h>

#include "generallensproblem.hh"
#include <dumux/common/start.hh>

////////////////////////
// the main function
////////////////////////
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0)
    {
        std::string errorMessageOut = "\nUsage: ";
        errorMessageOut += progName;
        errorMessageOut += " [options]\n";
        errorMessageOut += errorMsg;
        errorMessageOut += "\n\nThe List of Mandatory arguments for this program is:\n"
                "\t-TEnd                          The end of the simulation. [s] \n"
                "\t-DtInitial                     The initial timestep size. [s] \n"
                "\t-Grid.Cells                    Number of cells in (x,y)-direction [-]\n"
                "\t-Grid.UpperRight               Coordinates of upper right grid corner [m]\n";
        errorMessageOut += "\n\nThe Optional command line argument:\n"
                "\t-ModelType                     Can be: box (2p box model, default), cc (2p cc model), sequential (2p impes model)\n";
        std::cout << errorMessageOut << "\n";
    }
}

int main(int argc, char** argv)
{
    Dune::ParameterTree paramTree;
    std::string s(Dumux::readOptions_(argc, argv, paramTree));
    if (s.empty()) // everything was read correctly
    {
        // default model type is box
        const std::string modelType(paramTree.get<std::string>("ModelType", "box"));
        if (modelType == "box")
        {
            typedef TTAG(BoxGeneralLensProblem) ProblemTypeTag;
            GET_PROP(ProblemTypeTag, ParameterTree)::runTimeParams()["ModelType"] = "box";
            return Dumux::start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "cc")
        {
            typedef TTAG(CCGeneralLensProblem) ProblemTypeTag;
            return Dumux::start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "sequential")
        {
            typedef TTAG(SequentialGeneralLensProblem) ProblemTypeTag;
            return Dumux::start<ProblemTypeTag>(argc, argv, usage);
        }
        else
        {
            Dumux::ParameterException e("Unknown ModelType: " + modelType);
            std::cerr << e << ". Abort!" << std::endl
                      << "ModelType can be: box (2p box model), cc (2p cc model), sequential (2p impes model)" << std::endl;
            exit(1);
        }
    }
    else
        DUNE_THROW(Dumux::ParameterException, "Unknown command line option " << s);
}
