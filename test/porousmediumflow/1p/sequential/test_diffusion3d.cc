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

#if HAVE_DUNE_ALUGRID || HAVE_UG

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/io/grid/gridmanager.hh>

#include "test_diffusionproblem3d.hh"
#include "resultevaluation3d.hh"

namespace Dumux
{

/*!
 * \ingroup Start
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file.
 *
 * \param   argc    The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
 * \param   argv    The 'argv' argument of the main function: array of pointers to the argument strings
 * \param   usage   Callback function for printing the usage message
 */
int start(int argc,
          char **argv,
          void (*usage)(const char *, const std::string &))
{
    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    using TypeTag = Properties::TTag::DiffusionTest;

    auto defaultParams = [] (Dune::ParameterTree& p) {GetProp<TypeTag, Properties::ModelDefaultParameters>::defaultParams(p);};
    Parameters::init(argc, argv, defaultParams, usage);

    ////////////////////////////////////////////////////////////
    // get some optional parameters
    ////////////////////////////////////////////////////////////
    const int numRefine = getParam<int>("Grid.Refinement", 0);

    auto outputName = getParam<std::string>("Problem.OutputName", "");
    if (outputName.size())
        outputName.insert(0, "_");

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    try { gridManager.init(); }
    catch (...) {
        std::string usageMessage = "\n\t -> Creation of the grid failed! <- \n\n";
        usageMessage += defaultUsageMessage(argv[0]);
        usage(argv[0], usageMessage);
        throw;
    }

    // print grid info
    auto& grid = gridManager.grid();
    Dune::gridinfo(grid);

    //////////////////////////////////////////////////////////////////////
    // run the simulation
    /////////////////////////////////////////////////////////////////////

    Dune::Timer timer;
    bool consecutiveNumbering = true;

    //////////////////////////////////////////////////////////////////////
    // finite volume TPFA test problem
    /////////////////////////////////////////////////////////////////////

    using FVTypeTag = Properties::TTag::FVTest;
    using FVProblem = GetPropType<FVTypeTag, Properties::Problem>;
    auto fvDefaultParams = [] (Dune::ParameterTree& p) {GetProp<FVTypeTag, Properties::ModelDefaultParameters>::defaultParams(p);};
    Dumux::Parameters::init(argc, argv, fvDefaultParams, usage);

    std::shared_ptr<FVProblem> fvProblem = std::make_shared<FVProblem>(grid);
    // set output name
    std::string fvOutput = "test_diffusion3d_fv" + outputName;
    if (numRefine > 0)
        fvOutput += "_numRefine" + std::to_string(numRefine);
    fvProblem->setName(fvOutput);

    timer.reset();

    fvProblem->init();
    fvProblem->calculateFVVelocity();

    auto fvTime = timer.elapsed();

    fvProblem->writeOutput();

    Dumux::ResultEvaluation fvResult;
    fvResult.evaluate(grid.leafGridView(), *fvProblem, consecutiveNumbering);

    //////////////////////////////////////////////////////////////////////
    // finite volume MPFA-L test problem
    /////////////////////////////////////////////////////////////////////

    using MPFALTypeTag = Properties::TTag::FVMPFAL3DTestTypeTag;
    using MPFALProblem = GetPropType<MPFALTypeTag, Properties::Problem>;
    auto mpfalDefaultParams = [] (Dune::ParameterTree& p) {GetProp<MPFALTypeTag, Properties::ModelDefaultParameters>::defaultParams(p);};
    Dumux::Parameters::init(argc, argv, mpfalDefaultParams, usage);

    std::shared_ptr<MPFALProblem> mpfaProblem = std::make_shared<MPFALProblem>(grid);
    // set output name
    std::string fvmpfaOutput = "test_diffusion3d_fvmpfal" + outputName;
    if (numRefine > 0)
        fvmpfaOutput += "_numRefine" + std::to_string(numRefine);
    mpfaProblem->setName(fvmpfaOutput);

    timer.reset();

    mpfaProblem->init();
    mpfaProblem->calculateFVVelocity();

    auto mpfaTime = timer.elapsed();

    mpfaProblem->writeOutput();

    Dumux::ResultEvaluation mpfaResult;
    mpfaResult.evaluate(grid.leafGridView(), *mpfaProblem, consecutiveNumbering);

    //////////////////////////////////////////////////////////////////////
    // mimetic finite difference test problem
    /////////////////////////////////////////////////////////////////////

    using MimeticTypeTag = Properties::TTag::MimeticTest;
    using MimeticProblem = GetPropType<MimeticTypeTag, Properties::Problem>;
    auto mimeticDefaultParams = [] (Dune::ParameterTree& p) {GetProp<MimeticTypeTag, Properties::ModelDefaultParameters>::defaultParams(p);};
    Dumux::Parameters::init(argc, argv, mimeticDefaultParams, usage);

    std::shared_ptr<MimeticProblem> mimeticProblem = std::make_shared<MimeticProblem>(grid);
    // set output name
    std::string mimeticOutput = "test_diffusion3d_mimetic" + outputName;
    if (numRefine > 0)
        mimeticOutput += "_numRefine" + std::to_string(numRefine);
    mimeticProblem->setName(mimeticOutput);

    timer.reset();

    mimeticProblem->init();
    mimeticProblem->calculateFVVelocity();

    auto mimeticTime = timer.elapsed();

    mimeticProblem->writeOutput();

    Dumux::ResultEvaluation mimeticResult;
    mimeticResult.evaluate(grid.leafGridView(), *mimeticProblem, consecutiveNumbering);

    //////////////////////////////////////////////////////////////////////
    // print results to command line
    /////////////////////////////////////////////////////////////////////

    std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
    std::cout.precision(2);
    std::cout << "\t\t pErrorL2 \t pInnerErrorL2 \t vErrorL2 \t vInnerErrorL2 \t hMax \t\t pMin\t\t pMax\t\t pMinExact\t pMaxExact\t time"
              << std::endl;
    std::cout << "2pfa\t\t " << fvResult.relativeL2Error << "\t " << fvResult.relativeL2ErrorIn << "\t "
              << fvResult.ervell2 << "\t " << fvResult.ervell2In << "\t " << fvResult.hMax << "\t " << fvResult.uMin
              << "\t " << fvResult.uMax << "\t " << fvResult.uMinExact << "\t " << fvResult.uMaxExact << "\t "
              << fvTime << std::endl;
    std::cout << "mpfa-l\t " << mpfaResult.relativeL2Error << "\t " << mpfaResult.relativeL2ErrorIn << "\t "
              << mpfaResult.ervell2 << "\t " << mpfaResult.ervell2In << "\t " << mpfaResult.hMax << "\t "
              << mpfaResult.uMin << "\t " << mpfaResult.uMax << "\t " << mpfaResult.uMinExact << "\t "
              << mpfaResult.uMaxExact << "\t " << mpfaTime << std::endl;
    std::cout << "mimetic\t " << mimeticResult.relativeL2Error << "\t " << mimeticResult.relativeL2ErrorIn << "\t "
              << mimeticResult.ervell2 << "\t " << mimeticResult.ervell2In << "\t " << mimeticResult.hMax << "\t "
              << mimeticResult.uMin << "\t " << mimeticResult.uMax << "\t " << mimeticResult.uMinExact << "\t "
              << mimeticResult.uMaxExact << "\t " << mimeticTime << std::endl;

    return 0;
}

} // end namespace Dumux

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
                    errorMessageOut += "\n\nThe list of mandatory arguments for this program is:\n"
                                        "\t-Grid.File                      Name of the file containing the grid \n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

int main(int argc, char** argv)
{
    return Dumux::start(argc, argv, usage);
}
#else
int main()
{
#warning You need to have dune-ALUGrid or UG installed to run this test
    std::cerr << "You need to have dune-ALUGrid or UG installed to run this test\n";
    return 77;
}
#endif // HAVE_DUNE_ALUGRID || HAVE_UG
