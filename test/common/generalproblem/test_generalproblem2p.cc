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
#include "config.h"

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
                "\t-Grid.NumberOfCellsX           Resolution in x-direction [-]\n"
                "\t-Grid.NumberOfCellsY           Resolution in y-direction [-]\n"
                "\t-Grid.UpperRightX              Dimension of the grid [m]\n"
                "\t-Grid.UpperRightY              Dimension of the grid [m]\n";
        errorMessageOut += "\n\nThe Optional command line argument:\n"
                "\t-ModelType                     Can be: Box (2p box model), Decoupled (2p impes model),\n";
        std::cout << errorMessageOut << "\n";
    }
}

int main(int argc, char** argv)
{
    Dune::ParameterTree paramTree;
    std::string s(Dumux::readOptions_(argc, argv, paramTree));
    if (s.empty())
    {
        if (paramTree.hasKey("ModelType"))
        {
            std::string modelType(paramTree.get<std::string>("ModelType"));

            if (modelType == "Box")
            {
                typedef TTAG(BoxGeneralLensProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                ParamTree::tree()["Problem.OutputfileName"] = "generallens_box";
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used box 2p model\n";
                return startReturn;
            }
            else if (modelType == "CC")
            {
                typedef TTAG(CCGeneralLensProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                ParamTree::tree()["Problem.OutputfileName"] = "generallens_cc";
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used cc 2p model\n";
                return startReturn;
            }
            else if (modelType == "Decoupled")
            {
                typedef TTAG(DecoupledGeneralLensProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                ParamTree::tree()["Problem.OutputfileName"] = "generallens_decoupled";
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used decoupled 2p model\n";
                return startReturn;
            }
            else
            {
                typedef TTAG(BoxGeneralLensProblem) ProblemTypeTag;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Unknwon model type "<<modelType<<" specified\n";
                std::cout<<"Default to box model\n";
                return startReturn;
            }
        }
        else
        {
            typedef TTAG(BoxGeneralLensProblem) ProblemTypeTag;
            int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
            std::cout<<"######################################################\n";
            std::cout<<"No model type specified\n";
            std::cout<<"Default to box model\n";
            return startReturn;
        }
    }
    else
    {
        typedef TTAG(BoxGeneralLensProblem) ProblemTypeTag;
        int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
        std::cout<<"######################################################\n";
        std::cout<<s<<" is not a valid model type specification!\n";
        std::cout<<"Default to box model\n";
        return startReturn;
    }
}

