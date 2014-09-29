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
#include "config.h"

#include <dumux/common/start.hh>
#include <dumux/io/interfacemeshcreator.hh>
#include "evaporationatmosphereproblem.hh"

/*!
 * \brief Print a usage string for simulations.
 *
 * \param progName The name of the executable
 */
void printUsage(const char *progName)
{
    std::cout << "usage: " << progName
            << " [--restart restartTime] -parameterFile evaporationatmosphere.input\n";
    exit(1);
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

        // read a --my-opt=VALUE option. This gets transformed
        // into the parameter "MyOpt" with the value being "VALUE"
        if (argv[i][1] == '-') {
            // the syntax --my-opt=VALUE is deprecated and will be removed after DuMuX 2.6
            std::cout << std::endl
                << "Warning: the syntax --my-opt=VALUE is deprecated and will be removed after DuMuX 2.6"
                << std::endl << std::endl;

            std::string s(argv[i] + 2);
            // There is nothing after the '='
            if (s.size() == 0 || s[0] == '=')
            {
                std::ostringstream oss;
                oss << "\n -> Parameter name of argument " << i << " (='" << argv[i] << "')"
                    << " is empty. <- \n\n\n\n";
                return oss.str();
            }

            // capitalize the first character
            s[0] = toupper(s[0]);

            // parse argument
            unsigned int j = 0;
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
                else if (s[j] == '.') {
                    // we encountered a '.' character, indicating
                    // the end of a group name, and need
                    // to captitalize the following character
                    if (s.size() == j)
                    {
                        std::ostringstream oss;
                        oss << "\n -> Parameter name of argument " << i << " ('" << argv[i] << "')"
                            << " is invalid (ends with a '.' character). <- \n\n\n\n";
                        return oss.str();
                    }
                    s[j+1] = toupper(s[j+1]);
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
            // read a -MyOpt VALUE option
            paramName = argv[i] + 1;

            if (argc == i + 1 || argv[i+1][0] == '-') {
                std::ostringstream oss;
                oss << "\n -> No argument given for parameter '" << argv[i] << "'! <- \n\n\n\n";
                return oss.str();
            }

            paramValue = argv[i+1];
            ++i; // In the case of '-MyOpt VALUE' each pair counts as two arguments
        }

        // Put the key=value pair into the parameter tree
        paramTree[paramName] = paramValue;
    }
    return "";
}


template <class TypeTag>
int start_(int argc,
           char **argv)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef Dune::GridPtr<Grid> GridPointer;

    ////////////////////////////////////////////////////////////
    // Load the input parameters
    ////////////////////////////////////////////////////////////

//    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
//    Dune::ParameterTreeParser::readOptions(argc, argv, ParameterTree::tree());

    // fill the parameter tree with the options from the command line
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    std::string s = readOptions_(argc, argv, ParameterTree::tree());

    if (ParameterTree::tree().hasKey("ParameterFile") or argc==1)
    {
        // read input file, but do not overwrite options specified
        // on the command line, since the latter have precedence.
        std::string inputFileName ;
        if(argc==1) // if there are no arguments given (and there is a file ./<programname>.input) we use it as input file
        {
            std::cout<< "\nNo parameter file given. \n"
                     << "Defaulting to '"
                     << argv[0]
                     << ".input' for input file.\n";
            inputFileName = argv[0];
            inputFileName += ".input";
        }
        else
            inputFileName = GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile); // otherwise we read from the command line

        std::ifstream parameterFile;

        // check whether the parameter file exists.
        parameterFile.open(inputFileName.c_str());
        if (not parameterFile.is_open()){
            std::cout<< "\n\t -> Could not open file"
                     << inputFileName
                     << ". <- \n\n\n\n";
            printUsage(argv[0]);
            return 1;
        }
        parameterFile.close();

        Dune::ParameterTreeParser::readINITree(inputFileName,
                                               ParameterTree::tree(),
                                               /*overwrite=*/false);
    }

    // initialize MPI, finalize is done automatically on exit
    static Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // define the problem dimensions
    const int dim=2;

    // deal with the restart stuff
    int argIdx = 1;
    bool restart = false;
    double tStart = 0.0;
    if (argc > 1 && std::string("--restart") == argv[argIdx])
    {
        restart = true;
        ++argIdx;

        std::istringstream(argv[argIdx++]) >> tStart;
    }

    std::string dgfFileName;
    Scalar dt, tEnd;
    Dune::FieldVector<int, dim> numElements;
    Scalar interfacePos, gradingFactor;
    int gridRefinement;
    bool useInterfaceMeshCreator;

    try
    {
        dgfFileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Grid, File);
        dt = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);
        tEnd = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, TEnd);
        numElements[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, CellsX);
        if (dim>1) numElements[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, CellsY);
        if (dim==3) numElements[2] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, CellsZ);
        interfacePos = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePos);
        gradingFactor = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, Grading);
        gridRefinement = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, Refinement);
        useInterfaceMeshCreator = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Grid, UseInterfaceMeshCreator);
    }
    catch (Dumux::ParameterException &e) {
        std::cerr << e << ". Abort!\n";
        exit(1) ;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        exit(1);
    }
    std::cout << "Starting with timestep size = " << dt << "s, simulation end = " << tEnd << "s\n";

    GridPointer gridPtr;
    if (useInterfaceMeshCreator)
    {
        Dumux::InterfaceMeshCreator<Grid> interfaceMeshCreator;
        gridPtr = interfaceMeshCreator.create(dgfFileName, numElements, interfacePos, gradingFactor);
    }
    else
        gridPtr = GridPointer(dgfFileName);

    if (gridRefinement)
    {
        Grid& grid =  *gridPtr;
        grid.globalRefine(gridRefinement);
    }

    if (mpiHelper.size() > 1) {
        if (!Dune::Capabilities::isParallel<Grid>::v) {
            std::cerr << "WARNING: THE PROGRAM IS STARTED USING MPI, BUT THE GRID IMPLEMENTATION\n"
                      << "         YOU HAVE CHOSEN IS NOT PARALLEL!\n";
        }
        (*gridPtr).loadBalance();
    }

    // Instantiate the time manager
    TimeManager timeManager;

    // instantiate coupled problem
    Problem problem(timeManager, gridPtr->leafGridView());
    Dumux::Parameters::print<TypeTag>();

    // run the simulation
    timeManager.init(problem,
                     tStart, // initial time
                     dt, // initial time step
                     tEnd, // final time
                     restart);

    // print all properties
    Dumux::Properties::print<TypeTag>();

    timeManager.run();

    return 0;
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
 */
template <class TypeTag>
int start(int argc, char **argv)
{
    try {
        return start_<TypeTag>(argc, argv);
    }
    catch (Dumux::ParameterException &e)
    {
       std::cerr << e << ". Abort!\n";
       printUsage(argv[0]);
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

int main(int argc, char** argv)
{
#if !HAVE_UG && !(HAVE_ALUGRID || HAVE_DUNE_ALUGRID)
    std::cout<<"Evaporation Atmosphere not built, needs either UG or ALU for the log mesh." << std::endl;
    return 77;
#else
    typedef TTAG(EvaporationAtmosphereProblem) ProblemTypeTag;
    return start<ProblemTypeTag>(argc, argv);//, usage);
#endif

}
