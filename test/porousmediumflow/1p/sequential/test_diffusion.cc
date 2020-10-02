// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup IMPETtests
 * \brief test for diffusion models
 */
#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager.hh>

#include "test_diffusionproblem.hh"
#include "resultevaluation.hh"

////////////////////////
// the main function
////////////////////////

void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " #refine [delta]\n";
                    errorMessageOut += errorMsg;

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::FVVelocity2PTestTypeTag;

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    auto defaultParams = [] (Dune::ParameterTree& p) {GetProp<TypeTag, Properties::ModelDefaultParameters>::defaultParams(p);};
    Dumux::Parameters::init(argc, argv, defaultParams, usage);

    ////////////////////////////////////////////////////////////
    // create the grid
    ////////////////////////////////////////////////////////////
    Dumux::GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    auto& grid = gridManager.grid();

    ////////////////////////////////////////////////////////////
    // instantiate and run the concrete problem
    ////////////////////////////////////////////////////////////
    Dune::Timer timer;
    bool consecutiveNumbering = true;

    using FVProblem = GetPropType<Properties::TTag::FVVelocity2PTestTypeTag, Properties::Problem>;
    FVProblem fvProblem(grid);
    fvProblem.setName("fvdiffusion");
    timer.reset();
    fvProblem.init();
    fvProblem.calculateFVVelocity();
    double fvTime = timer.elapsed();
    fvProblem.writeOutput();
    Dumux::ResultEvaluation fvResult;
    fvResult.evaluate(grid.leafGridView(), fvProblem, consecutiveNumbering);

    using MPFAOProblem = GetPropType<Properties::TTag::FVMPFAOVelocity2PTestTypeTag, Properties::Problem>;
    MPFAOProblem mpfaProblem(grid);
    mpfaProblem.setName("fvmpfaodiffusion");
    timer.reset();
    mpfaProblem.init();
    double mpfaTime = timer.elapsed();
    mpfaProblem.writeOutput();
    Dumux::ResultEvaluation mpfaResult;
    mpfaResult.evaluate(grid.leafGridView(), mpfaProblem, consecutiveNumbering);

    using MimeticProblem = GetPropType<Properties::TTag::MimeticPressure2PTestTypeTag, Properties::Problem>;
    MimeticProblem mimeticProblem(grid);
    mimeticProblem.setName("mimeticdiffusion");
    timer.reset();
    mimeticProblem.init();
    double mimeticTime = timer.elapsed();
    mimeticProblem.writeOutput();
    Dumux::ResultEvaluation mimeticResult;
    mimeticResult.evaluate(grid.leafGridView(), mimeticProblem, consecutiveNumbering);

    std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
    std::cout.precision(2);
    std::cout << "\t error press \t error grad\t sumflux\t erflm\t\t uMin\t\t uMax\t\t time" << std::endl;
    std::cout << "2pfa\t " << fvResult.relativeL2Error << "\t " << fvResult.ergrad << "\t " << fvResult.sumflux
                    << "\t " << fvResult.erflm << "\t " << fvResult.uMin
                    << "\t " << fvResult.uMax << "\t " << fvTime << std::endl;
    std::cout << "mpfa-o\t " << mpfaResult.relativeL2Error << "\t " << mpfaResult.ergrad
                    << "\t " << mpfaResult.sumflux << "\t " << mpfaResult.erflm
                    << "\t " << mpfaResult.uMin << "\t " << mpfaResult.uMax << "\t " << mpfaTime << std::endl;
    std::cout << "mimetic\t " << mimeticResult.relativeL2Error << "\t " << mimeticResult.ergrad
                    << "\t " << mimeticResult.sumflux << "\t " << mimeticResult.erflm
                    << "\t " << mimeticResult.uMin << "\t " << mimeticResult.uMax << "\t " << mimeticTime << std::endl;



    return 0;
}
