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
 * \ingroup Common
 * \brief Provides a few default main functions for convenience.
 */
#ifndef DUMUX_COMMON_START_HH
#define DUMUX_COMMON_START_HH

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/io/grid/gridmanager.hh>

#warning "start.hh is deprecated. Use new style main files see e.g. /test/porousmediumflow/1p."

namespace Dumux {

/*!
 * \ingroup Common
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param   argc    The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
 * \param   argv    The 'argv' argument of the main function: array of pointers to the argument strings
 * \param   usage   Callback function for printing the usage message
 */
template <class TypeTag>
int start_(int argc,
           char **argv,
           void (*usage)(const char *, const std::string &))
{
    // some aliases for better readability
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    auto defaultParams = [] (Dune::ParameterTree& p) {GetProp<TypeTag, Properties::ModelDefaultParameters>::defaultParams(p);};
    Parameters::init(argc, argv, defaultParams, usage);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    //////////////////////////////////////////////////////////////////////
    // run the simulation
    /////////////////////////////////////////////////////////////////////

    // read the initial time step and the end time (mandatory parameters)
    auto tEnd = getParam<Scalar>("TimeManager.TEnd");
    auto dt = getParam<Scalar>("TimeManager.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    bool restart = false;
    Scalar restartTime = 0;
    if (hasParam("Restart") || hasParam("TimeManager.Restart"))
    {
        restart = true;
        restartTime =  getParam<Scalar>("TimeManager.Restart");
    }

    // instantiate and run the problem
    TimeManager timeManager;
    Problem problem(timeManager, gridManager.grid());
    timeManager.init(problem, restartTime, dt, tEnd, restart);
    timeManager.run();

    // print dumux end message and maybe the parameters for debugging
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}

/*!
 * \ingroup Common
 *
 * \brief Provides a main function with error handling
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc  The number of command line arguments of the program
 * \param argv  The contents of the command line arguments of the program
 * \param usage Callback function for printing the usage message
 */
template <class TypeTag>
int start(int argc,
          char **argv,
          void (*usage)(const char *, const std::string &))
{
    try {
        return start_<TypeTag>(argc, argv, usage);
    }
    catch (ParameterException &e) {
        std::cerr << std::endl << e << ". Abort!" << std::endl;
        return 1;
    }
    catch (Dune::DGFException & e) {
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << std::endl;
    return 2;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 3;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        return 4;
    }
}

} // end namespace Dumux

#endif
