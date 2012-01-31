// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Provides a few default main functions for convenience.
 */
#ifndef DUMUX_START_HH
#define DUMUX_START_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/common/mpihelper.hh>
#include <iostream>

#include <dune/common/version.hh>
#include <dune/common/parametertreeparser.hh>

#include <sys/ptrace.h>

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(TimeManager);
}

/*!
 * \brief Print a usage string for simulations using
 *        Dumux::startFromDGF() as their main() function.
 *
 * \param progname The name of the executable
 */
void printUsageDGF(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] gridfile.dgf tEnd dt\n";
    exit(1);
}

/*!
 * \brief Print a usage string for simulations using
 *        Dumux::startWithGrid() as their main() function.
 *
 * \param progname The name of the executable
 */
void printUsageGrid(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] tEnd dt\n";
    exit(1);
}

/*!
 * \brief Print a usage string for simulations using
 *        Dumux::startWithGrid() as their main() function.
 *
 * \param progname The name of the executable
 */
void printUsageInputFile(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] inputfile\n";
    exit(1);
}

/*!
 * \ingroup Start
 * \brief Provides a default main function for simulations requiring
 *        only a single DGF file as their grid specification.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc  The 'argc' argument of the main function
 * \param argv  The 'argv' argument of the main function
 */
template <class TypeTag>
int startFromDGF(int argc, char **argv)
{
#ifdef NDEBUG
    try {
#endif

        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
        typedef Dune::GridPtr<Grid> GridPointer;

        // initialize MPI, finalize is done automatically on exit
        static Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc < 4)
            printUsageDGF(argv[0]);

        // deal with the restart stuff
        int argIdx = 1;
        bool restart = false;
        double startTime = 0;
        if (std::string("--restart") == argv[argIdx]) {
            restart = true;
            ++argIdx;

            std::istringstream(argv[argIdx++]) >> startTime;
        }

        if (argc - argIdx != 3) {
            printUsageDGF(argv[0]);
        }

        double tEnd, dt;
        const char *dgfFileName = argv[argIdx++];
        std::istringstream(argv[argIdx++]) >> tEnd;
        std::istringstream(argv[argIdx++]) >> dt;

        // create grid
        // -> load the grid from file
        GridPointer gridPtr(dgfFileName);
        if (mpiHelper.size() > 1) {
            if (!Dune::Capabilities::isParallel<Grid>::v) {
                std::cerr << "DUMUX WARNING: THE PROGRAM IS STARTED USING MPI, BUT THE GRID IMPLEMENTATION\n"
                          << "               YOU HAVE CHOSEN IS NOT PARALLEL!\n";
            }
            gridPtr.loadBalance();
        }

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager, gridPtr->leafView());
        timeManager.init(problem, startTime, dt, tEnd, restart);

        // print all properties
        Dumux::Properties::print<TypeTag>();

        timeManager.run();
        return 0;

#ifdef NDEBUG
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
    }
#endif

    return 3;
}


/*!
 * \ingroup Start
 * \brief Provides a default main function for simulations which
 *        create the grid themselves but do not require any other
 *        parameters.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param grid  The grid used by the simulation
 * \param argc  The 'argc' argument of the main function
 * \param argv  The 'argv' argument of the main function
 */
template <class TypeTag>
int startWithGrid(const typename GET_PROP_TYPE(TypeTag, Grid) &grid,
                  int argc,
                  char **argv)
{
#ifdef NDEBUG
    try {
#endif

        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

        // parse the command line arguments for the program
        if (argc < 3)
            printUsageGrid(argv[0]);

        // deal with the restart stuff
        int argIdx = 1;
        bool restart = false;
        double startTime = 0;
        if (std::string("--restart") == argv[argIdx]) {
            restart = true;
            ++argIdx;

            std::istringstream(argv[argIdx++]) >> startTime;
        }

        if (argc - argIdx != 2) {
            printUsageGrid(argv[0]);
        }

        double tEnd, dt;
        std::istringstream(argv[argIdx++]) >> tEnd;
        std::istringstream(argv[argIdx++]) >> dt;

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager, grid.leafView());
        timeManager.init(problem, startTime, dt, tEnd, restart);

        // print all properties
        Dumux::Properties::print<TypeTag>();

        timeManager.run();
        return 0;

#ifdef NDEBUG
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
    }
#endif

    return 3;
}

std::string readOptions_(int argc, char **argv, Dune::ParameterTree &paramTree)
{
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            std::ostringstream oss;
            oss << "Command line argument " << i << " (='" << argv[i] << "') is invalid.";
            return oss.str();
        }
        
        std::string paramName, paramValue;

        if (argv[i][1] == '-') {
            // read a --my-opt=abc option. This gets transformed
            // into the parameter "MyOpt" with the value being
            // "abc"
            std::string s(argv[i] + 2);
            if (s.size() == 0 || s[0] == '=')
            {
                std::ostringstream oss;
                oss << "Parameter name of argument " << i << " (='" << argv[i] << "')"
                    << " is empty.";
                return oss.str();
            }

            
            // capitalize first letter
            s[0] = std::toupper(s[0]);

            // parse argument
            int j = 0;
            while (true) {
                if (j >= s.size()) {
                    // encountered the end of the string, i.e. we
                    // have a parameter where the argument is empty
                    paramName = s;
                    paramValue = "";
                    break;
                }
                else if (s[j] == '=') {
                    // we encountered a '=' character. everything
                    // before is the name of the parameter,
                    // everything after is the value.
                    paramName = s.substr(0, j);
                    paramValue = s.substr(j+1);
                    break;
                }
                else if (s[j] == '-') {
                    // remove all "-" characters and capitalize the
                    // character after them
                    s.erase(j, 1);
                    if (s.size() == j)
                    {
                        std::ostringstream oss;
                        oss << "Parameter name of argument " << i << " ('" << argv[i] << "')"
                            << " is invalid (ends with a '-' character).";
                        return oss.str();
                    }
                    else if (s[j] == '-')
                    {
                        std::ostringstream oss;
                        oss << "Malformed parameter name name in argument " << i << " ('" << argv[i] << "'): "
                            << "'--' in parameter name.";
                        return oss.str();
                    }
                    s[j] = toupper(s[j]);
                };

                ++j;
            }
        }
        else {
            // read a -MyOpt abc option
            paramName = argv[i] + 2;
            
            if (argc == i + 1 || argv[i+1][0] == '-') {
                std::ostringstream oss;
                oss << "No argument given for parameter '" << argv[i] << "'!";
                return oss.str();
            }
            
            paramValue = argv[i+1];
            ++i;
        }

        // Put the key=value pair into the parameter tree
        paramTree[paramName] = paramValue;
    }
    return "";
}

template <class TypeTag>
int startWithParameters_(int argc,
                         char **argv, 
                         void (*usage)(const char *, const std::string &))
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
        
    // initialize MPI, finalize is done automatically on exit
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
    
    ////////////////////////////////////////////////////////////
    // parse the command line arguments
    ////////////////////////////////////////////////////////////

    // check whether the user did not specify any parameter. in this
    // case print the usage message
    if (argc == 1) {
        usage(argv[0], "");
        return 0;
    }


    // check whether the user wanted to see the help message
    for (int i = 1; i < argc; ++i) {
        if (std::string("--help") == argv[i] || std::string("-h") == argv[i])
        {
            usage(argv[0], "");
            return 0;
        }
    }
    
    // fill the parameter tree with the options from the command line
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;
    std::string s = readOptions_(argc, argv, Params::tree());
    if (!s.empty()) {
        usage(argv[0], s);
        return 1;
    };
    
    if (Params::tree().hasKey("OptsFile")) {
        // read input file, but do not overwrite options specified
        // on the command line, since the latter have preceedence.
        std::string inputFileName = GET_RUNTIME_PARAM(TypeTag, std::string, OptsFile);
        Dune::ParameterTreeParser::readINITree(inputFileName, 
                                               Params::tree(), 
                                               /*overwrite=*/false);
    }


    bool printProps = true;
    if (Params::tree().hasKey("PrintProperties"))
        printProps = GET_RUNTIME_PARAM(TypeTag, bool, PrintProperies);

    if (printProps && mpiHelper.rank() == 0) {
        Dumux::Properties::print<TypeTag>();
    }
    
    // deal with the restart stuff
    bool restart = false;
    Scalar restartTime = 0;
    if (Params::tree().hasKey("Restart")) {
        restart = true;
        restartTime = GET_RUNTIME_PARAM(TypeTag, Scalar, Restart);
    }
    
    // read the PrintParams parameter
    bool printParams = true;
    if (Params::tree().hasKey("PrintParameters"))
        printParams = GET_RUNTIME_PARAM(TypeTag, bool, PrintParameters);
    
    try { GridCreator::makeGrid(); }
    catch (...) { usage(argv[1], "Creation of the grid failed!"); throw; }
    
    // read the initial time step and the end time
    double tEnd;
    double dt;

    try { tEnd = GET_RUNTIME_PARAM(TypeTag, Scalar, TEnd); }
    catch (...) { usage(argv[1], "Mandatory parameter '--t-end' not specified!"); throw; }

    try { dt = GET_RUNTIME_PARAM(TypeTag, Scalar, DtInitial); }
    catch (...) { usage(argv[1], "Mandatory parameter '--dt-initial' not specified!"); throw; }

    // instantiate and run the concrete problem
    TimeManager timeManager;
    Problem problem(timeManager);
    timeManager.init(problem, 0, dt, tEnd, restart);
    if (restart)
        problem.restart(restartTime);
    timeManager.run();

    if (printParams && mpiHelper.rank() == 0) {
        Dumux::Parameters::print<TypeTag>();
    }
        
    return 0;
}

/*!
 * \ingroup Start
 *
 * \brief Returns true if and only if a debugger is attached to the simulation.
 */
bool inDebugger()
{
    // valgrind seems to have a problem with ptrace, so we behave as
    // if no debugger is present in this case...
    if (Valgrind::Running())
        return false;

    return ptrace(PTRACE_TRACEME, 0, NULL, 0) == -1;
}


/*!
 * \ingroup Start
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc The number of command line arguments of the program
 * \param argv The contents of the command line arguments of the program
 * \param usage Callback function for printing the usage message
 */
template <class TypeTag>
int startWithParameters(int argc, 
                        char **argv,
                        void (*usage)(const char *, const std::string &))
{
    if (!inDebugger()) {
        try {
            return startWithParameters_<TypeTag>(argc, argv, usage);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            return 1;
        }
        catch (Dune::Exception &e) {
            std::cerr << "Dune reported error: " << e << std::endl;
            return 2;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            return 3;
        }
    }
    else
        return startWithParameters_<TypeTag>(argc, argv, usage);
}

} // namespace Dumux

#endif
