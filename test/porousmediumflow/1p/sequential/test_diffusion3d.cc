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

#include <dumux/common/start.hh>

#include "test_diffusionproblem3d.hh"
#include "resultevaluation3d.hh"

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
    try {
        typedef TTAG(DiffusionTestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef GET_PROP(TypeTag, ParameterTree) ParameterTree;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // fill the parameter tree with the options from the command line
        std::string s = Dumux::readOptions_(argc, argv, ParameterTree::tree());
        if (!s.empty()) {
            usage(argv[0], s);
            return 1;
        }

        // obtain the name of the parameter file
        std::string parameterFileName;
        if (ParameterTree::tree().hasKey("ParameterFile"))
        {
            // set the name to the one provided by the user
            parameterFileName = GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile);
        }
        else // otherwise we read from the command line
        {
            // set the name to the default ./<programname>.input
            parameterFileName = argv[0];
            parameterFileName += ".input";
        }

        Dune::ParameterTreeParser::readINITree(parameterFileName,
                                               ParameterTree::tree(),
                                               /*overwrite=*/false);

        int numRefine = 0;
        if (ParameterTree::tree().hasKey("Grid.RefinementRatio"))
        numRefine = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, RefinementRatio);

        std::string outputName("");
        if (ParameterTree::tree().hasKey("Problem.OutputName"))
        {
            outputName += "_";
            outputName += GET_RUNTIME_PARAM(TypeTag, std::string, OutputName);
        }

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        std::string gridFileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Grid, File);
        Dune::GridPtr < Grid > grid(gridFileName);
        grid->globalRefine(numRefine);

        Dune::gridinfo (*grid);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        Dune::Timer timer;
        bool consecutiveNumbering = true;

        typedef TTAG (FVTestProblem)
        FVTypeTag;
        typedef GET_PROP_TYPE(FVTypeTag, Problem) FVProblem;
        typedef GET_PROP(FVTypeTag, ParameterTree) FVParameterTree;
        Dune::ParameterTreeParser::readINITree(parameterFileName, FVParameterTree::tree());
        FVProblem *fvProblem = new FVProblem(grid->leafGridView());

        std::string fvOutput("test_diffusion3d_fv");
        fvOutput += outputName;
        if (numRefine > 0)
        {
            char refine[128];
            sprintf(refine, "_numRefine%d", numRefine);
            fvOutput += refine;
        }
        fvProblem->setName(fvOutput.c_str());
        timer.reset();
        fvProblem->init();
        fvProblem->calculateFVVelocity();
        double fvTime = timer.elapsed();
        fvProblem->writeOutput();
        Dumux::ResultEvaluation fvResult;
        fvResult.evaluate(grid->leafGridView(), *fvProblem, consecutiveNumbering);
        delete fvProblem;

        typedef TTAG (FVMPFAL3DTestProblem)
        MPFALTypeTag;
        typedef GET_PROP_TYPE(MPFALTypeTag, Problem) MPFALProblem;
        typedef GET_PROP(MPFALTypeTag, ParameterTree) MPFALParameterTree;
        Dune::ParameterTreeParser::readINITree(parameterFileName, MPFALParameterTree::tree());
        MPFALProblem *mpfaProblem = new MPFALProblem(grid->leafGridView());

        std::string fvmpfaOutput("test_diffusion3d_fvmpfal");
        fvmpfaOutput += outputName;
        if (numRefine > 0)
        {
            char refine[128];
            sprintf(refine, "_numRefine%d", numRefine);
            fvmpfaOutput += refine;
        }
        mpfaProblem->setName(fvmpfaOutput.c_str());
        timer.reset();
        mpfaProblem->init();
        double mpfaTime = timer.elapsed();
        mpfaProblem->writeOutput();
        Dumux::ResultEvaluation mpfaResult;
        mpfaResult.evaluate(grid->leafGridView(), *mpfaProblem, consecutiveNumbering);
        delete mpfaProblem;

        typedef TTAG (MimeticTestProblem)
        MimeticTypeTag;
        typedef GET_PROP_TYPE(MimeticTypeTag, Problem) MimeticProblem;
        typedef GET_PROP(MimeticTypeTag, ParameterTree) MimeticParameterTree;
        Dune::ParameterTreeParser::readINITree(parameterFileName, MimeticParameterTree::tree());
        MimeticProblem *mimeticProblem = new MimeticProblem(grid->leafGridView());

        std::string mimeticOutput("test_diffusion3d_mimetic");
        mimeticOutput += outputName;
        if (numRefine > 0)
        {
            char refine[128];
            sprintf(refine, "_numRefine%d", numRefine);
            mimeticOutput += refine;
        }
        mimeticProblem->setName(mimeticOutput.c_str());
        timer.reset();
        mimeticProblem->init();
        double mimeticTime = timer.elapsed();
        mimeticProblem->writeOutput();
        Dumux::ResultEvaluation mimeticResult;
        mimeticResult.evaluate(grid->leafGridView(), *mimeticProblem, consecutiveNumbering);
        delete mimeticProblem;

        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        std::cout.precision(2);
        std::cout
                << "\t pErrorL2 \t pInnerErrorL2 \t vErrorL2 \t vInnerErrorL2 \t hMax \t\t pMin\t\t pMax\t\t pMinExact\t pMaxExact\t time"
                << std::endl;
        std::cout << "2pfa\t " << fvResult.relativeL2Error << "\t " << fvResult.relativeL2ErrorIn << "\t "
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
    } catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...)
    {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

    return 3;
}
#else
int main()
{
#warning You need to have dune-ALUGrid or UG installed to run this test
    std::cerr << "You need to have dune-ALUGrid or UG installed to run this test\n";
    return 77;
}
#endif // HAVE_DUNE_ALUGRID || HAVE_UG
