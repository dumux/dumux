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
 * \brief Provides a few default main functions for convenience.
 */
#ifndef DUMUX_START_HH
#define DUMUX_START_HH

#include <ctime>
#include <iostream>

#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/dgfparser.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/dgf.hh>
#endif // HAVE_DUNE_ALUGRID

#include "propertysystem.hh"
#include "parameters.hh"
#include "valgrind.hh"

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(TimeManager);
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


/*!
 * \ingroup Start
 *
 * \brief Provides a general text block, that is part of error/ help messages.
 *
 * \return The string that is the help / error message.
 */
std::string usageTextBlock()
{
    return  "Options usually are parameters given to the simulation, \n"
            "and have to be specified with this syntax: \n"
            "\t-GroupName.ParameterName VALUE, for example -TimeManager.TEnd 100\n"
            "\n"
            "If -ParameterFile is specified, parameters can also be defined there. In this case,\n"
            "lines of the form \n"
            "GroupName.ParameterName = VALUE # comment \n"
            "have to be used. More conveniently, group names can be specified in square brackets, \n"
            "such that each following parameter name belongs to that group, \n"
            "[GroupName] \n"
            "ParameterName = VALUE \n"
            "See test/boxmodels/2p/test_2p.input for an example of a parameter file \n"
            "and the Dune documentation of ParameterTreeParser for the format specification. \n"
            "\n"
            "Parameters specified on the command line have priority over those in the parameter file.\n"
            "If no parameter file name is given, './<programname>.input' is chosen as default.\n"
            "\n"
            "Important options include:\n"
            "\t-h, --help                        Print this usage message and exit\n"
            "\t-PrintParameters [true|false]     Print the run-time modifiable parameters _after_ \n"
            "\t                                  the simulation [default: true]\n"
            "\t-PrintProperties [true|false]     Print the compile-time parameters _before_ \n"
            "\t                                  the simulation [default: false]\n"
            "\t-ParameterFile FILENAME           File with parameter definitions\n"
            "\t-TimeManager.Restart RESTARTTIME  Restart simulation from a restart file\n"
            "\n"
            "\n";
}

/*!
 * \ingroup Start
 * \brief Dumux start and end message.
 *
 * Selects random messages to write out at the start and end of a simulation run.
 * \param start Indicates if the simulation is starting (=true) or ending (=false).
 */
void dumuxMessage_(bool start)
{
    static int dice_;
    if(start)
    {
        std::srand(std::time(0));
        // roll the dice to decide which start message will be displayed:
        dice_ = std::rand() % 12; // has to be n+1
    }


    std::cout << "" << std::endl;
    switch (dice_)
    {
        case 0:
            if(start)
                std::cout << "Welcome aboard DuMuX airlines. Please fasten your seatbelts! "
                          << "Emergency exits are near the time integration." << std::endl;
            else
                std::cout << "We hope that you enjoyed simulating with us " << std::endl
                          << "and that you will choose us next time, too." << std::endl;
        break;
        case 1:
            if(start)
                std::cout << "Let's get the cow off the ice." << std::endl;
            else
                std::cout << "DuMuX got the cow off the ice." << std::endl;
        break;
        case 2:
            if(start)
                std::cout << "Science, my lad, is made up of mistakes, but they are "
                          << "mistakes which it is useful to make, because they lead little "
                          << "by little to the truth." << std::endl
                          << " - Jules Verne, A journey to the center of the earth" << std::endl;
            else
                std::cout << "[We see that] science is eminently perfectible, and that each theory has "
                          << "constantly to give way to a fresh one." << std::endl
                          << " - Jules Verne, Journey to the Center of the Earth" << std::endl;

        break;
        case 3:
            if(start)
                std::cout << "Wherever he saw a hole he always wanted to know the depth of it. "
                          << "To him this was important." << std::endl
                          << " - Jules Verne, A journey to the center of the earth" << std::endl;
            else
                std::cout << "We may brave human laws, but we cannot resist natural ones." << std::endl
                          << " - Jules Verne, 20,000 Leagues Under the Sea" << std::endl;
        break;
        case 4:
            if(start)
                std::cout << "Silence - to delight Bernd." << std::endl;
            else
                std::cout << "\n " << std::endl;
        break;
        case 5:
                std::cout << "Don't panic... !" << std::endl;
        break;
        case 6:
            if(start)
                std::cout << "You idiot! You signed the order to destroy Earth!" << std::endl
                          << " - Douglas Adams, HGttG" << std::endl;
            else
                std::cout << "Marvin: I've been talking to the main computer.\n Arthur: And?" << std::endl
                          << "Marvin: It hates me." << std::endl
                          << " - Douglas Adams, HGttG" << std::endl;
        break;
        case 7:
            if(start)
                std::cout << "In the beginning the Universe was created. This has made a lot of "
                          << "people very angry and has been widely regarded as a bad move.!" << std::endl
                          << " - Douglas Adams, HGttG " << std::endl;
            else
                std::cout << "Forty-two. I checked it very thoroughly, and that quite definitely is the answer. I think "
                          << "the problem, to be quite honest with you, is that you\'ve never actually known what the question is." << std::endl
                          << " - Douglas Adams, HGttG " << std::endl;
        break;
        case 8:
            std::cout << "                ##                  @@@@          @   @         @ @" << std::endl;
            std::cout << "             ###   #                @   @         @@ @@          @ " << std::endl;
            std::cout << "           ##       #               @   @  @   @  @ @ @  @   @  @ @" << std::endl;
            std::cout << "         ##          #              @   @  @   @  @   @  @   @     " << std::endl;
            std::cout << "        #             #             @@@@    @@@   @   @   @@@   2.6" << std::endl;
            std::cout << "       #               #                                           " << std::endl;
            std::cout << "      #                 #                                          " << std::endl;
            std::cout << "     #                   ##        %%%                             " << std::endl;
            std::cout << "    #                      ###    %   %  %%     %%                 " << std::endl;
            std::cout << "####                          #%%%     %%  %%%%%  %%%%%%%%%%%%%%%%%" << std::endl;
        break;
        case 9:
            std::cout << "###         #   #        # #                            " << std::endl;
            std::cout << "#  #  #  #  ## ##  #  #   #                             " << std::endl;
            std::cout << "#  #  #  #  # # #  #  #  # #                            " << std::endl;
            std::cout << "###    ##   #   #   ##                                  " << std::endl;
            std::cout << "                                                        " << std::endl;
            std::cout << "Dune for Multi-{ Phase,                                 " << std::endl;
            std::cout << "                 Component,                             " << std::endl;
            std::cout << "                 Scale,                                 " << std::endl;
            std::cout << "                 Physics,                               " << std::endl;
            std::cout << "                 ...} flow and transport in porous media" << std::endl;
        break;
        case 10:
          if(start)
            {
                std::cout << "Elliot Carver: Mr. Jones, are we ready to release our new software?" << std::endl;
                std::cout << "Jones: Yes, sir. As requested, it's full of bugs, which means people will be forced to upgrade for years." << std::endl;
                std::cout << " - James Bond, Tomorrow Never Dies" << std::endl;
            }
            else
            {
                std::cout << "Elliot Carver: Outstanding." << std::endl;
                std::cout << " - James Bond, Tomorrow Never Dies" << std::endl;
            }
        break;
        case 11:
          if(start)
            {
                std::cout << "Chuck Norris has successfully compiled DuMuX." << std::endl;
            }
            else
            {
                std::cout << "Chuck Norris has compiled DuMuX even two times in row!" << std::endl;
            }
        break;
        //TODO: If you add a case, you have to increase the modulus at the beginning of the function!

        default:    // silence to delight Bernd
            return ;
    }
    std::cout << std::endl;
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
int start_(int argc,
           char **argv,
           void (*usage)(const char *, const std::string &))
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // Set by default (dumux/common/basicproperties.hh) to DgfGridCreator (dumux/io/dgfgridcreator.hh)
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    // initialize MPI, finalize is done automatically on exit
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        dumuxMessage_(true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments
    ////////////////////////////////////////////////////////////

    // check whether the user did not specify any parameter. in this
    // case print the usage message
    if (argc == 1) {
        if (mpiHelper.size() > 0)
            std::cout << "Rank " << mpiHelper.rank() << ": ";
        std::cout << "No parameter file given. "
                  << "Defaulting to '"
                  << argv[0]
                  << ".input' for input file.\n";
        std::ifstream parameterFile;
        // check whether the parameter file exists.
        std::string defaultName = argv[0];
        defaultName += ".input";
        parameterFile.open(defaultName.c_str());
        if (not parameterFile.is_open()){
            if (mpiHelper.size() > 0)
                std::cout << "Rank " << mpiHelper.rank() << ": ";
            std::cout << " -> Could not open file '"
                      << defaultName
                      << "'. <- \n\n\n\n";
            usage(argv[0], usageTextBlock());
            return 1;
        }
        parameterFile.close();
    }


    // check whether the user wanted to see the help message
    if (mpiHelper.rank() == 0)
    {
        for (int i = 1; i < argc; ++i)
        {
            if (std::string("--help") == argv[i] || std::string("-h") == argv[i])
            {
                usage(argv[0], usageTextBlock());
                return 0;
            }
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

    // obtain the name of the parameter file
    std::string parameterFileName;
    if (ParameterTree::tree().hasKey("ParameterFile"))
    {
        // set the name to the one provided by the user
        parameterFileName = GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile); // otherwise we read from the command line        
    }
    else 
    {
        // set the name to the default ./<programname>.input
        parameterFileName = argv[0];
        parameterFileName += ".input";
    }
    
    // open and check whether the parameter file exists.
    std::ifstream parameterFile(parameterFileName.c_str());
    if (not parameterFile.is_open()) {
        // if the name of the file has been specified, 
        // this must be an error. Otherwise proceed.
        if (ParameterTree::tree().hasKey("ParameterFile"))
        {
            std::cout << "\n\t -> Could not open file '"
                      << parameterFileName
                      << "'. <- \n\n\n\n";
            usage(argv[0], usageTextBlock());
            return 1;
        }
    }
    else 
    {
        // read parameters from the file without overwriting
        Dune::ParameterTreeParser::readINITree(parameterFileName,
                                               ParameterTree::tree(),
                                               /*overwrite=*/false);
    }
    parameterFile.close();

    bool printProps = false;
    if (ParameterTree::tree().hasKey("PrintProperties") 
        || ParameterTree::tree().hasKey("TimeManager.PrintProperties"))
        printProps = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, TimeManager, PrintProperties);

    if (printProps && mpiHelper.rank() == 0) {
        Dumux::Properties::print<TypeTag>();
    }

    // deal with the restart stuff
    bool restart = false;
    Scalar restartTime = 0;
    if (ParameterTree::tree().hasKey("Restart") 
        || ParameterTree::tree().hasKey("TimeManager.Restart")) {
        restart = true;
        restartTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, Restart);
    }

    // read the PrintParams parameter
    bool printParams = true;
    if (ParameterTree::tree().hasKey("PrintParameters") 
        || ParameterTree::tree().hasKey("TimeManager.PrintParameters"))
        printParams = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, TimeManager, PrintParameters);

    // try to create a grid (from the given grid file)
    try { GridCreator::makeGrid(); }
    catch (...) {
        std::string usageMessage = "\n\t -> Creation of the grid failed! <- \n\n\n\n";
        usageMessage += usageTextBlock();
        usage(argv[0], usageMessage);
        throw;
    }
    GridCreator::loadBalance();

    // read the initial time step and the end time
    double tEnd;
    double dt;
 
    try { tEnd = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, TEnd); }
    catch (...) {
        std::string usageMessage = "\n\t -> Mandatory parameter 'TimeManager.TEnd'"
                                   " not specified! <- \n\n\n\n";
        usageMessage += usageTextBlock();
        usage(argv[0], usageMessage);
        throw;
    }

    try { dt = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial); }
    catch (...) {
        std::string usageMessage = "\n\t -> Mandatory parameter 'TimeManager.DtInitial'"
                                   " not specified! <- \n\n\n\n";
        usageMessage += usageTextBlock();
        usage(argv[0], usageMessage);
        throw;
    }

    // instantiate and run the concrete problem
    TimeManager timeManager;
    Problem problem(timeManager, GridCreator::grid().leafGridView());
    timeManager.init(problem, restartTime, dt, tEnd, restart);
    timeManager.run();
    // print dumux end message
    if (mpiHelper.rank() == 0)
        dumuxMessage_(false);

    if (printParams && mpiHelper.rank() == 0) {
        Dumux::Parameters::print<TypeTag>();
    }
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
    catch (Dumux::ParameterException &e) {
        std::cerr << e << ". Abort!\n";
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

} // namespace Dumux

#endif
