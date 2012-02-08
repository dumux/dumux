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

#include <iostream>
#include <sys/ptrace.h>

#include "propertysystem.hh"
#include "parameters.hh"
#include "valgrind.hh"

#include <dune/common/mpihelper.hh>
#include <dune/common/version.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/dgfparser.hh>

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Scalar);
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

/*!
 * \ingroup Start
 * \brief Read the command line arguments and write them into the parameter tree.
 *        Do some syntax checks.
 *
 * \param   argc      The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
 * \param   argv      The 'argv' argument of the main function: array of pointers to the argument strings
 * \param   paramTree The parameterTree. It can be filled from an input file or the command line.
 * \return            Empty string if everything worked out. Otherwise the thing that could not be read.
 */
std::string readOptions_(int argc, char **argv, Dune::ParameterTree &paramTree)
{
    // All command line options need to start with '-'
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            std::ostringstream oss;
            oss << "\n -> Command line argument " << i << " (='" << argv[i] << "') is invalid. <- \n\n\n\n";
            return oss.str();
        }

        std::string paramName, paramValue;

        // read a --my-opt=abc option. This gets transformed
        // into the parameter "myOpt" with the value being
        // "abc"
        if (argv[i][1] == '-') {
            std::string s(argv[i] + 2);
            // There is nothing after the '='
            if (s.size() == 0 || s[0] == '=')
            {
                std::ostringstream oss;
                oss << "\n -> Parameter name of argument " << i << " (='" << argv[i] << "')"
                    << " is empty. <- \n\n\n\n";
                return oss.str();
            }

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
                        oss << "\n -> Parameter name of argument " << i << " ('" << argv[i] << "')"
                            << " is invalid (ends with a '-' character). <- \n\n\n\n";
                        return oss.str();
                    }
                    else if (s[j] == '-')
                    {
                        std::ostringstream oss;
                        oss << "\n -> Malformed parameter name name in argument " << i << " ('" << argv[i] << "'): "
                            << "'--' in parameter name. <- \n\n\n\n";
                        return oss.str();
                    }
                    s[j] = toupper(s[j]);
                }

                ++j;
            }
        }
        else {
            // read a -myOpt abc option
            paramName = argv[i] + 1;

            if (argc == i + 1 || argv[i+1][0] == '-') {
                std::ostringstream oss;
                oss << "\n -> No argument given for parameter '" << argv[i] << "'! <- \n\n\n\n";
                return oss.str();
            }

            paramValue = argv[i+1];
            ++i; // In the case of '-myOpt abc' each pair counts as two arguments
        }

        // Put the key=value pair into the parameter tree
        paramTree[paramName] = paramValue;
    }
    return "";
}


/*!
 * \ingroup Start
 *
 * \brief Provides a general text block, that is part of error/ help messages.
 *
 * \return The string that is the help / error message.
 */
std::string usageTextBlock()
{
    return  "Options have to be specified with this syntax: \n"
            "\t-tEnd ENDTIME                    The time of the end of the simulation [s]\n"
            "Alternativ supported syntax:\n"
            "\t--t-end=ENDTIME                  The time of the end of the simulation [s]\n"
            "\n"
            "If -parameterFile is specified parameters can also be defined there. In this case,\n"
            "camel case is used for the parameters (e.g.: tEnd=100). \n"
            "\n"
            "Parameters specified on the command line have priority over those in the parameter file.\n"
            "\n"
            "Important optional options include:\n"
            "\t-h, --help                        Print this usage message and exit\n"
            "\t-PrintParameters [true|false]     Print the run-time modifiable parameters _after_ \n"
            "\t                                  the simulation [default: true]\n"
            "\t-PrintProperties [true|false]     Print the compile-time parameters _before_ \n"
            "\t                                  the simulation [default: false]\n"
            "\t-parameterFile FILENAME           File with parameter definitions\n"
            "\t-restart RESTARTTIME              Restart simulation from a restart file\n"
            "\n"
            "For the case of no arguments given, the input parameter file is expected to be named './parameter.input' \n"
            "\n";
}

/*!
 * \ingroup Start
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
int startWithParameters_(int argc,
                         char **argv,
                         void (*usage)(const char *, const std::string &))
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator; // Set by default (dumux/common/basicproperties.hh) to DgfGridCreator (dumux/common/dgfgridcreator.hh)
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
        std::cout<< "\nNo parameter file given. \n"
                 << "Defaulting to './parameter.input' for input file.\n";
        std::ifstream parameterFile;
        // check whether the parameter file exists.
        parameterFile.open("parameter.input");
        if (not parameterFile.is_open()){
            std::cout<< "\n\t -> Could not open file './parameter.input'. <- \n\n\n\n";
            usage(argv[0], usageTextBlock());
            return 1;
        }
        parameterFile.close();
    }


    // check whether the user wanted to see the help message
    for (int i = 1; i < argc; ++i) {
        if (std::string("--help") == argv[i] || std::string("-h") == argv[i])
        {
            usage(argv[0], usageTextBlock());
            return 0;
        }
    }

    // fill the parameter tree with the options from the command line
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    std::string s = readOptions_(argc, argv, ParameterTree::tree());
    if (!s.empty()) {
        std::string usageMessage = s ;
                    usageMessage += usageTextBlock();
        usage(argv[0], usageMessage);
        return 1;
    }

    if (ParameterTree::tree().hasKey("parameterFile") or argc==1) {
        // read input file, but do not overwrite options specified
        // on the command line, since the latter have precedence.
        std::string inputFileName ;
        if(argc==1) // if there are no arguments given (and there is a file ./parameter.input) we use it as input file
            inputFileName="parameter.input";
        else
            inputFileName = GET_RUNTIME_PARAM(TypeTag, std::string, parameterFile); // otherwise we read from the command line

        std::ifstream parameterFile;

        // check whether the parameter file exists.
        parameterFile.open(inputFileName.c_str());
        if (not parameterFile.is_open()){
            std::cout<< "\n\t -> Could not open file"
                     << inputFileName
                     << ". <- \n\n\n\n";
            usage(argv[0], usageTextBlock());
            return 1;
        }
        parameterFile.close();

        Dune::ParameterTreeParser::readINITree(inputFileName,
                                               ParameterTree::tree(),
                                               /*overwrite=*/false);
    }

    bool printProps = false;
    if (ParameterTree::tree().hasKey("PrintProperties"))
        printProps = GET_RUNTIME_PARAM(TypeTag, bool, PrintProperties);

    if (printProps && mpiHelper.rank() == 0) {
        Dumux::Properties::print<TypeTag>();
    }

    // deal with the restart stuff
    bool restart = false;
    Scalar restartTime = 0;
    if (ParameterTree::tree().hasKey("restart")) {
        restart = true;
        restartTime = GET_RUNTIME_PARAM(TypeTag, Scalar, restart);
    }

    // read the PrintParams parameter
    bool printParams = true;
    if (ParameterTree::tree().hasKey("PrintParameters"))
        printParams = GET_RUNTIME_PARAM(TypeTag, bool, PrintParameters);

    // try to create a grid (from the given grid file)
    try { GridCreator::makeGrid(); }
    catch (...) {
        std::string usageMessage = "\n\t -> Creation of the grid failed! <- \n\n\n\n";
                    usageMessage += usageTextBlock();
        usage(argv[0], usageMessage);
        throw;
    }

    // read the initial time step and the end time
    double tEnd;
    double dt;

    try { tEnd = GET_RUNTIME_PARAM(TypeTag, Scalar, tEnd); }
    catch (...) {
        std::string usageMessage = "\n\t -> Mandatory parameter '--t-end' not specified! <- \n\n\n\n";
                    usageMessage += usageTextBlock();
        usage(argv[0], usageMessage);
        throw;
    }

    try { dt = GET_RUNTIME_PARAM(TypeTag, Scalar, dtInitial); }
    catch (...) {
        std::string usageMessage = "\n\t -> Mandatory parameter '--dt-initial' not specified! <- \n\n\n\n";
                    usageMessage += usageTextBlock();
        usage(argv[0], usageMessage);
        throw;
    }

    // instantiate and run the concrete problem
    TimeManager timeManager;
    Problem problem(timeManager, GridCreator::grid().leafView());
    timeManager.init(problem, restartTime, dt, tEnd, restart);
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
 *        In this function only the differentiation between debugger
 *        or not is made.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc  The number of command line arguments of the program
 * \param argv  The contents of the command line arguments of the program
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
