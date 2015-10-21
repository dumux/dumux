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
 * \brief Test for the sequential 2p models in 3d
 */
#define PROBLEM 2 // 0 = Buckley-Leverett, 1 = McWhorter, 2 = Nine-Spot

#include "config.h"

#if HAVE_DUNE_ALUGRID

#include "test_3d2pproblem.hh"
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
#if STRUCTUREDGRID
                "\t-Grid.NumberOfCellsX           Resolution in x-direction [-]\n"
                "\t-Grid.NumberOfCellsY           Resolution in y-direction [-]\n"
                "\t-Grid.NumberOfCellsZ           Resolution in y-direction [-]\n"
                "\t-Grid.UpperRightX              Dimension of the grid [m]\n"
                "\t-Grid.UpperRightY              Dimension of the grid [m]\n"
                "\t-Grid.UpperRightZ              Dimension of the grid [m]\n";
#else
                "\t-Grid.File                     Name of the grid file (*.dgf)";
#endif
        errorMessageOut += "\n\nThe Optional command line argument:\n"
                "\t-ModelType                     Can be: FV (standard finite volume), FVAdaptive (adaptive finite volume),\n"
#if PROBLEM == 1
                "\t                     MPFAL (MPFA l-method), MPFALAdaptive (adaptive MPFA l-method)\n";
#else
                "\t                     MPFAL (MPFA l-method), MPFALAdaptive (adaptive MPFA l-method)\n"
                "\t                     Mimetic (mimetic FDM), MimeticAdaptive (adaptive MFDM)\n";
#endif
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

            if (modelType == "FV")
            {
                typedef TTAG(FVTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used standard finite volume model\n";
                return startReturn;
            }
            else if (modelType == "FVAdaptive")
            {
                typedef TTAG(FVAdaptiveTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used adaptive finite volume model\n";
                return startReturn;
            }
            else if (modelType == "MPFAL")
            {
                typedef TTAG(MPFALTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used finite volume MPFA l-method model\n";
                return startReturn;
            }
            else if (modelType == "MPFALAdaptive")
            {
                typedef TTAG(MPFALAdaptiveTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used adaptive finite volume MPFA l-method model\n";
                return startReturn;
            }
#if PROBLEM != 1
            else if (modelType == "Mimetic")
            {
                typedef TTAG(MimeticTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used mimetic finite difference model\n";
                return startReturn;
            }
            else if (modelType == "MimeticAdaptive")
            {
                typedef TTAG(MimeticAdaptiveTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::runTimeParams();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Used adaptive mimetic finite difference model\n";
                return startReturn;
            }
#endif
            else
            {
                typedef TTAG(MPFALTwoPTestProblem) ProblemTypeTag;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
                std::cout<<"######################################################\n";
                std::cout<<"Unknwon model type "<<modelType<<" specified\n";
                std::cout<<"Default to finite volume MPFA l-method model\n";
                return startReturn;
            }
        }
        else
        {
            typedef TTAG(MPFALTwoPTestProblem) ProblemTypeTag;
            int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
            std::cout<<"######################################################\n";
            std::cout<<"No model type specified\n";
            std::cout<<"Default to finite volume MPFA l-method model\n";
            return startReturn;
        }
    }
    else
    {
        typedef TTAG(MPFALTwoPTestProblem) ProblemTypeTag;
        int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv, usage);
        std::cout<<"######################################################\n";
        std::cout<<s<<" is not a valid model type specification!\n";
        std::cout<<"Default to finite volume MPFA l-method model\n";
        return startReturn;
    }

    return 0;
}

#else

#include <iostream>

int main()
{
#warning You need to have dune-ALUGrid installed to run this test
    std::cerr << "You need to have dune-ALUGrid installed to run this test\n";
    return 77;
}
#endif
