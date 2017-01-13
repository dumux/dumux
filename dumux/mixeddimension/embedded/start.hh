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
 * \ingroup EmbeddedCoupling
 * \brief Provides a starting algorithm (a default main function) for multidimension problems
 */
#ifndef DUMUX_START_EMBEDDED_MIXEDDIMENSION_HH
#define DUMUX_START_EMBEDDED_MIXEDDIMENSION_HH

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/parameterparser.hh>

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(TimeManager);
NEW_PROP_TAG(LowDimProblemTypeTag);
NEW_PROP_TAG(BulkProblemTypeTag);
}

template <class TypeTag>
int start_(int argc,
           char **argv,
           void (*usage)(const char *, const std::string &))
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ParameterTree = typename GET_PROP(TypeTag, ParameterTree);
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);
    using BulkGridCreator = typename GET_PROP_TYPE(BulkProblemTypeTag, GridCreator);
    using LowDimGridCreator = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridCreator);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // if the user just wanted to see the help / usage message show usage and stop program
    if(!ParameterParser::parseCommandLineArguments(argc, argv, ParameterTree::tree(), usage))
    {
        usage(argv[0], defaultUsageMessage(argv[0]));
        return 0;
    }

    // parse the input file into the parameter tree
    // check first if the user provided an input file through the command line, if not use the default
    const auto parameterFileName = ParameterTree::tree().hasKey("ParameterFile") ? GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile) : "";
    ParameterParser::parseInputFile(argc, argv, ParameterTree::tree(), parameterFileName, usage);

    ////////////////////////////////////////////////////////////
    // check for some user debugging parameters
    ////////////////////////////////////////////////////////////

    bool printProps = false; // per default don't print all properties
    if (ParameterTree::tree().hasKey("PrintProperties") || ParameterTree::tree().hasKey("TimeManager.PrintProperties"))
        printProps = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, TimeManager, PrintProperties);

    if (printProps && mpiHelper.rank() == 0)
        Properties::print<TypeTag>();

    bool printParams = true; // per default print all properties
    if (ParameterTree::tree().hasKey("PrintParameters") || ParameterTree::tree().hasKey("TimeManager.PrintParameters"))
        printParams = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, TimeManager, PrintParameters);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    try { BulkGridCreator::makeGrid(); }
    catch (...) {
        std::string usageMessage = "\n\t -> Creation of the bulk grid failed! <- \n\n";
        usageMessage += defaultUsageMessage(argv[0]);
        usage(argv[0], usageMessage);
        throw;
    }
    BulkGridCreator::loadBalance();

    try { LowDimGridCreator::makeGrid(); }
    catch (...) {
        std::string usageMessage = "\n\t -> Creation of the low dim grid failed! <- \n\n";
        usageMessage += defaultUsageMessage(argv[0]);
        usage(argv[0], usageMessage);
        throw;
    }
    LowDimGridCreator::loadBalance();


    //////////////////////////////////////////////////////////////////////
    // run the simulation
    /////////////////////////////////////////////////////////////////////

    // read the initial time step and the end time (mandatory parameters)
    auto tEnd = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, TEnd);
    auto dt = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);

    // check if we are about to restart a previously interrupted simulation
    bool restart = false;
    Scalar restartTime = 0;
    if (ParameterTree::tree().hasKey("Restart") || ParameterTree::tree().hasKey("TimeManager.Restart"))
    {
        restart = true;
        restartTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, Restart);
    }

    // instantiate and run the concrete problem
    TimeManager timeManager;
    Problem problem(timeManager, BulkGridCreator::grid().leafGridView(), LowDimGridCreator::grid().leafGridView());
    timeManager.init(problem, restartTime, dt, tEnd, restart);
    timeManager.run();

    // print dumux end message and maybe the parameters for debugging
    if (mpiHelper.rank() == 0)
    {
        DumuxMessage::print(/*firstCall=*/false);

        if (printParams)
            Parameters::print<TypeTag>();
    }

    return 0;
}

/*!
 * \ingroup Start
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
        Parameters::print<TypeTag>();
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
}

} // end namespace Dumux

#endif
