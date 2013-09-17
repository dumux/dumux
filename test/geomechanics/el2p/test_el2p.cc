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
 * \brief test for the el2p CC model
 */
#include "config.h"

#if HAVE_DUNE_PDELAB

// Check if DUNE-PDELab has been patched for our needs.
#ifdef DUNE_PDELAB_IS_PATCHED_FOR_DUMUX

#include "el2pproblem.hh"
#include <dune/common/precision.hh>
#include <dumux/common/start.hh>

int main(int argc, char** argv) {
    try {
        typedef TTAG(El2P_TestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::GridPtr<Grid> GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        if (argc < 5) {
            std::cout<<"usage: "<<argv[0]<<" grid tEnd_Initialization tEnd dt "<<std::endl;
            return 1;
        }

        int argPos = 1;
        double tEnd(0), dt(0), tInitEnd(0);

        // get the grid name
        const char *dgfFileName = argv[argPos++];
        // get the end of the initialization period
        std::istringstream(argv[argPos++]) >> tInitEnd;
        // get the end of the real simulation period
        std::istringstream(argv[argPos++]) >> tEnd;
        // get the initial time step (applied for initialization and for real simulation)
        std::istringstream(argv[argPos++]) >> dt;

        // load the grid from file
        GridPointer gridPtr(dgfFileName);
        (*gridPtr).loadBalance();

        // Instantiate the time manager
        TimeManager timeManager;

        // instantiate problem
        Problem problem(timeManager, gridPtr->leafView(),tInitEnd);

        // set the initial approximated hydrostatic pressure distribution
        // based on an averaged brine density
        // or based on a pressure polynomial
        problem.initializePressure();

        // start initialization run to initialize the pressure field correctly
        timeManager.init(problem, 0.0, // initial time
                dt, // initial time step
                tInitEnd); // end of initialization period
        std::cout<<"tInit: "<<tInitEnd<<" tEnd: "<<tEnd<<" dt: "<<dt<<std::endl;
        timeManager.run();

        // for the real simulation the coupling between mass balances and momentum equation
        // is turned on
        problem.setCoupled(true);
        // pressure field resulting from the initialization period is applied for the initial
        // and the Dirichlet boundary conditions
        problem.setPressure();
        // output is written
        problem.setOutput(true);
        // run the real simulation
        timeManager.init(problem, tInitEnd, // initial time
                dt, // initial time step
                tEnd + tInitEnd); // final time

        timeManager.run();

        return 0;
    } catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}

#else // DUNE_PDELAB_IS_PATCHED_FOR_DUMUX

#warning You need to have  a patched dune-pdelab to run this test, see ../../../patches/README for details.

#include <iostream>

int main()
{
    std::cerr << "You need to have a patched dune-pdelab to run this test, "
                 "see ../../../patches/README for details." << std::endl;;
    return 77;
}

#endif // DUNE_PDELAB_IS_PATCHED_FOR_DUMUX

#else // HAVE_DUNE_PDELAB

#warning You need to have dune-pdelab installed and patched to run this test.

#include <iostream>

int main()
{
    std::cerr << "You need to have dune-pdelab installed and patched to run this test.\n";
    return 77;
}
#endif // HAVE_DUNE_PDELAB
