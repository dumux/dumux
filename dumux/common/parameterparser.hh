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
 * \brief Provides a class parsing parameters from command line and from input files
 */
#ifndef DUMUX_PARAMETER_PARSER_HH
#define DUMUX_PARAMETER_PARSER_HH

#include <iostream>

#include <dune/common/parametertreeparser.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/exceptions.hh>

namespace Dumux
{

/*!
 * \ingroup Start
 * \brief Parses parameters in the command line and input files
 */
template<class TypeTag>
class ParameterParser
{

public:

    /*!
     * \brief Read the command line arguments and write them into the parameter tree.
     *        Do some syntax checks. Check also if an input file was specified in the command line.
     *        In parallel the parameter are currently read on all processors.
     *
     * \param   argc      The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
     * \param   argv      The 'argv' argument of the main function: array of pointers to the argument strings
     * \param   params    A parameter tree. It can be filled from an input file or the command line.
     * \param   usage     Callback function for printing the usage message
     * \return  bool      True if everything succeeded, false if the help message has to be shown
     */
    static bool parseCommandLineArguments(int argc,
                                          char **argv,
                                          Dune::ParameterTree &params,
                                          void (*usage)(const char *, const std::string &) = [](const char *, const std::string &){})
    {
        const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

        // check whether the user did not specify any parameter then we are done
        if (argc == 1)
            return true;

        // check whether the user wanted to see the help message
        if (mpiHelper.rank() == 0)
            for (int i = 1; i < argc; ++i)
                if (std::string("--help") == argv[i] || std::string("-h") == argv[i])
                    return false;

        // fill the parameter tree with the parameters from the command line
        readOptions_(argc, argv, params);

        return true;
    }

    /*!
     * \brief Parse the input file. If the user didn't specify anything in the parameter tree (that contains
     *        command line options if parseCommandLineArguments was called before) we default to the
     *        program name + input. When calling this function we consider it an error if no input file is found.
     *
     * \param   argc      The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
     * \param   argv      The 'argv' argument of the main function: array of pointers to the argument strings
     * \param   params    A parameter tree. It can be filled from an input file or the command line.
     * \param   usage     Callback function for printing the usage message
     */
    static void parseInputFile(int argc,
                               char **argv,
                               Dune::ParameterTree &params,
                               void (*usage)(const char *, const std::string &) = [](const char *, const std::string &){})
    {
        const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

        std::string parameterFileName = "";

        // check the parameter tree for a user specified input file
        if (params.hasKey("ParameterFile"))
            // this is the only reason why this class needs a TypeTag -- sigh
            // if the parameter tree is used directly this appears in the unused properties list
            parameterFileName = GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile);

        // otherwise use the default
        else
        {
            if (mpiHelper.size() > 1)
                std::cout << "Rank " << mpiHelper.rank() << ": ";
            std::cout << "No parameter file given. "
                      << "Defaulting to '"
                      << argv[0]
                      << ".input' for input file.\n";

            parameterFileName = std::string(argv[0]) + ".input";
        }

        // open and check whether the parameter file exists.
        std::ifstream parameterFile(parameterFileName.c_str());
        if (!parameterFile.is_open())
        {
            if (mpiHelper.size() > 1)
                std::cout << "Rank " << mpiHelper.rank() << ": ";
            std::cout << " -> Could not open file '"
                      << parameterFileName
                      << "'. <- \n\n";

            usage(argv[0], defaultUsageMessage(argv[0]));

            DUNE_THROW(ParameterException, "Error opening input file " << parameterFileName << ".");
        }
        else
        {
            // read parameters from the file without overwriting
            Dune::ParameterTreeParser::readINITree(parameterFileName,
                                                   params,
                                                   /*overwrite=*/false);
        }
        parameterFile.close();
    }

private:
    /*!
     * \brief Read the command line arguments and write them into the parameter tree.
     *        Do some syntax checks.
     *
     * \param   argc      The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
     * \param   argv      The 'argv' argument of the main function: array of pointers to the argument strings
     * \param   paramTree The parameterTree. It can be filled from an input file or the command line.
     * \return            Empty string if everything worked out. Otherwise the thing that could not be read.
     */
    static void readOptions_(int argc, char **argv, Dune::ParameterTree &paramTree)
    {
        // All command line options need to start with '-'
        for (int i = 1; i < argc; ++i)
        {
            if (argv[i][0] != '-' && i == 1)
            {
                // try to pass first argument as parameter file
                paramTree["ParameterFile"] = argv[1];
                continue;
            }

            if (argv[i][0] != '-')
                DUNE_THROW(ParameterException, "-> Command line argument " << i << " (='" << argv[i] << "') is invalid. <-");

            if (i+1 == argc)
                DUNE_THROW(ParameterException, "-> No argument given for parameter '" << argv[i] << "'! <-");

            // read a -MyOpt VALUE option
            std::string paramName = argv[i] + 1;
            std::string paramValue = argv[i+1];
            ++i; // In the case of '-MyOpt VALUE' each pair counts as two arguments

            // Put the key=value pair into the parameter tree
            paramTree[paramName] = paramValue;
        }
    }
};

} // end namespace Dumux

#endif
