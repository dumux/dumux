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
 * \brief Test for the sequential 2p models
 */
#define PROBLEM 2 // 0 = Buckley-Leverett, 1 = McWhorter, 2 = 2D Lense problem

#include <config.h>

#if HAVE_UG

#include <dumux/common/start.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include "test_mpfa2pproblem.hh"

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
                "\t-Grid.Cells           Resolution in x- and y-direction [-]\n"
                "\t-Grid.UpperRight              Dimension of the grid [m]\n";
        errorMessageOut += "\n\nThe Optional command line argument:\n"
                "\t-ModelType                     Can be: FV (standard finite volume), FVAdaptive (adaptive finite volume),\n"
                "\t                     MPFAO (MPFA o-method), MPFAL (MPFA l-method), MPFALAdaptive (adaptive MPFA l-method)\n";
        std::cout << errorMessageOut << "\n";
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    try {
        auto defaultParams = [] (Dune::ParameterTree& p) {GET_PROP(TTAG(MPFATwoPTestTypeTag), ModelDefaultParameters)::defaultParams(p);};
        Parameters::init(argc, argv, defaultParams, usage);

        const auto modelType = getParam<std::string>("ModelType", "MPFAL");

        if (modelType == "FV")
        {
            using ProblemTypeTag = TTAG(FVTwoPTestTypeTag);
            std::cout<<"##########################################" << std::endl;
            std::cout<<"Standard finite volume TPFA model" << std::endl;
            std::cout<<"##########################################" << std::endl;
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "FVAdaptive")
        {
            using ProblemTypeTag = TTAG(FVAdaptiveTwoPTestTypeTag);
            std::cout<<"##########################################" << std::endl;
            std::cout<<"Adaptive finite volume TPFA model" << std::endl;
            std::cout<<"##########################################" << std::endl;
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "MPFAO")
        {
            using ProblemTypeTag = TTAG(MPFAOTwoPTestTypeTag);
            std::cout<<"##########################################" << std::endl;
            std::cout<<"Standard finite volume MPFA-O model" << std::endl;
            std::cout<<"##########################################" << std::endl;
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "MPFAL")
        {
            using ProblemTypeTag = TTAG(MPFALTwoPTestTypeTag);
            std::cout<<"##########################################" << std::endl;
            std::cout<<"Unknown model type " << modelType << ", default to" << std::endl;
            std::cout<<"Standard finite volume MPFA-L model" << std::endl;
            std::cout<<"##########################################" << std::endl;
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else if (modelType == "MPFALAdaptive")
        {
            using ProblemTypeTag = TTAG(MPFALAdaptiveTwoPTestTypeTag);
            std::cout<<"##########################################" << std::endl;
            std::cout<<"Adaptive finite volume MPFA-L model" << std::endl;
            std::cout<<"##########################################" << std::endl;
            return start<ProblemTypeTag>(argc, argv, usage);
        }
        else
        {
            using ProblemTypeTag = TTAG(MPFAOTwoPTestTypeTag);
            std::cout<<"##########################################" << std::endl;
            std::cout<<"Unknown model type " << modelType << ", default to" << std::endl;
            std::cout<<"Standard finite volume MPFA-O model" << std::endl;
            std::cout<<"##########################################" << std::endl;
            return start<ProblemTypeTag>(argc, argv, usage);
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

    return 0;
}

#else

#include <iostream>

int main()
{
#warning You need to have UGGrid to run this test
    std::cerr << "You need to have UGGrid to run this test\n";
    return 77;
}
#endif
