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

#include <dumux/common/start.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/parameterparser.hh>

#include "generallensproblem.hh"

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
    using namespace Dumux;

    try {

        Dune::ParameterTree paramTree;
        // if the user just wanted to see the help / usage message show usage and stop program
        if(!ParameterParser::parseCommandLineArguments(argc, argv, paramTree, usage))
        {
            usage(argv[0], defaultUsageMessage(argv[0]));
            return 0;
        }

        // default model type is box
        const std::string modelType(paramTree.get<std::string>("ModelType", "box"));
        if (modelType == "box")
        {
            using ProblemTypeTag = TTAG(BoxGeneralLensProblem);
            // avoid unused parameter message
            GET_PROP(ProblemTypeTag, ParameterTree)::runTimeParams()["ModelType"] = "box";
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "cc")
        {
            using ProblemTypeTag = TTAG(CCGeneralLensProblem);
            // avoid unused parameter message
            GET_PROP(ProblemTypeTag, ParameterTree)::runTimeParams()["ModelType"] = "cc";
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "sequential")
        {
            using ProblemTypeTag = TTAG(SequentialGeneralLensProblem);
            // avoid unused parameter message
            GET_PROP(ProblemTypeTag, ParameterTree)::runTimeParams()["ModelType"] = "sequential";
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else
        {
            std::cerr << ParameterException("Unknown ModelType: " + modelType) << ". Abort!" << std::endl
                      << "ModelType can be: box (2p box model), cc (2p cc model), sequential (2p impes model)" << std::endl;
            return 1;
        }
    }
    catch (ParameterException &e) {
        std::cerr << std::endl << e << ". Abort!" << std::endl;
        return 1;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 3;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        return 4;
    }

} // end main
