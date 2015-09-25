// $Id$
/*****************************************************************************
 *   Copyright (C) 20010 by Markus Wolff                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \ingroup IMPETtests
 * \brief test for diffusion models
 */
#include "config.h"
#include <iostream>
#include <boost/format.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/gridinfo.hh>

#include "test_diffusion_problem.hh"
#include "benchmarkresult.hh"

////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << boost::format("usage: %s #refine [delta]\n")%progname;
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(DiffusionTestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        static const int dim = Grid::dimension;
        typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc != 2 && argc != 3)
            usage(argv[0]);

        int numRefine;
        std::istringstream(argv[1]) >> numRefine;

        double delta = 1e-3;
        if (argc == 3)
            std::istringstream(argv[2]) >> delta;

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        Dune::FieldVector<int,dim> N(1);
        GlobalPosition L(0.0);
        GlobalPosition H(1.0);
        Grid grid(N,L,H);
        grid.globalRefine(numRefine);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        Dune::Timer timer;
        bool consecutiveNumbering = true;

        typedef GET_PROP_TYPE(TTAG(FVVelocity2PTestProblem), PTAG(Problem)) FVProblem;
        FVProblem fvProblem(grid.leafView(), delta);
        timer.reset();
        fvProblem.init();
        fvProblem.model().calculateVelocity();
        double fvTime = timer.elapsed();
        fvProblem.writeOutput();
        Dumux::ResultEvaluation fvResult;
        fvResult.evaluate(grid.leafView(), fvProblem, fvProblem.variables().pressure(), fvProblem.variables().velocity(), consecutiveNumbering);

        typedef GET_PROP_TYPE(TTAG(FVMPFAOVelocity2PTestProblem), PTAG(Problem)) MPFAOProblem;
        MPFAOProblem mpfaProblem(grid.leafView(), delta);
        timer.reset();
        mpfaProblem.init();
        mpfaProblem.model().calculateVelocity();
        double mpfaTime = timer.elapsed();
        mpfaProblem.writeOutput();
        Dumux::ResultEvaluation mpfaResult;
        mpfaResult.evaluate(grid.leafView(), mpfaProblem, mpfaProblem.variables().pressure(), mpfaProblem.variables().velocity(), consecutiveNumbering);

        typedef GET_PROP_TYPE(TTAG(MimeticPressure2PTestProblem), PTAG(Problem)) MimeticProblem;
        MimeticProblem mimeticProblem(grid.leafView(), delta);
        timer.reset();
        mimeticProblem.init();
        mimeticProblem.model().calculateVelocity();
        double mimeticTime = timer.elapsed();
        mimeticProblem.writeOutput();
        Dumux::ResultEvaluation mimeticResult;
        mimeticResult.evaluate(grid.leafView(), mimeticProblem, mimeticProblem.variables().pressure(), mimeticProblem.variables().velocity(), consecutiveNumbering);

        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        std::cout.precision(2);
        std::cout << "\t error press \t error grad\t sumflux\t erflm\t\t uMin\t\t uMax\t\t time" << std::endl;
        std::cout << "2pfa\t " << fvResult.relativeL2Error << "\t " << fvResult.ergrad << "\t " << fvResult.sumflux
                        << "\t " << fvResult.erflm << "\t " << fvResult.uMin << "\t " << fvResult.uMax << "\t " << fvTime << std::endl;
        std::cout << "mpfa-o\t " << mpfaResult.relativeL2Error << "\t " << mpfaResult.ergrad << "\t " << mpfaResult.sumflux
                        << "\t " << mpfaResult.erflm << "\t " << mpfaResult.uMin << "\t " << mpfaResult.uMax << "\t " << mpfaTime << std::endl;
        std::cout << "mimetic\t " << mimeticResult.relativeL2Error << "\t " << mimeticResult.ergrad << "\t " << mimeticResult.sumflux
                        << "\t " << mimeticResult.erflm << "\t " << mimeticResult.uMin << "\t " << mimeticResult.uMax << "\t " << mimeticTime << std::endl;



        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

    return 3;
}
