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
 * \ingroup Parameter
 * \file
 *
 * \brief The infrastructure to retrieve run-time parameters from
 *        Dune::ParameterTrees with the defaul value taken from the
 *        property system.
 */
#ifndef DUMUX_PARAMETERS_HH
#define DUMUX_PARAMETERS_HH

#include <iostream>
#include <list>
#include <sstream>
#include <unordered_map>
#include <fstream>

#include <dune/common/parametertree.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/loggingparametertree.hh>

/*!
 * \ingroup Parameter
 * \brief Retrieve a runtime parameter which _does_ have a default value taken from
 *        the Dumux property system.
 *
 * Example:
 *
 * \code
 * // -> retrieves scalar value UpwindWeight, default
 * // is taken from the property UpwindWeight
 * GET_PARAM(TypeTag, Scalar, UpwindWeight);
 * \endcode
 */
// #define GET_PARAM(TypeTag, ParamType, ParamName)
//     ::Dumux::template getParam_UsingDeprecatedMacro<ParamType>(std::string(#ParamName), GET_PROP_VALUE(TypeTag, ParamName))

/*!
 * \ingroup Parameter
 * \brief Retrieve a runtime parameter which _does_ have a default value taken from
 *        the Dumux property system.
 *
 * The third argument is group name which must be the prefix to the
 * property name which provides the default value for the parameter
 *
 * Example:
 *
 * \code
 * // -> retrieves Boolean value Newton.WriteConvergence, default
 * // is taken from the property NewtonWriteConvergence
 * GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence);
 * \endcode
 */
#define GET_PARAM_FROM_GROUP(TypeTag, ParamType, GroupName, ParamName)  \
    ::Dumux::template getParam_UsingDeprecatedMacro<ParamType>(std::string(#GroupName) + "." + std::string(#ParamName), GET_PROP_VALUE(TypeTag, GroupName##ParamName))


/*!
 * \ingroup Parameter
 * \brief Retrieve a runtime parameter which _does not_ have a default value taken from
 *        the Dumux property system.
 *
 * Example:
 *
 * \code
 * // -> retrieves global integer value NumberOfCellsX
 * GET_RUNTIME_PARAM(TypeTag, int, NumberOfCellsX);
 * \endcode
 */
#define GET_RUNTIME_PARAM(TypeTag, ParamType, ParamName) \
    ::Dumux::template getParam_UsingDeprecatedMacro<ParamType>(std::string(#ParamName))

/*!
 * \ingroup Parameter
 * \brief Retrieve a runtime parameter which _does not_ have a default value taken from
 *        the Dumux property system.
 *
 * The third argument is the complete parameter name, which has to be a c-string.
 * This allows to use string variables as parameter name.
 *
 * Example with a temporary c-string:
 *
 * \code
 * // -> retrieves global integer value "NumberOfCellsX" which is
 * // located in the parameter group "Grid"
 * GET_RUNTIME_PARAM_CSTRING(TypeTag, int, "Grid.NumberOfCellsX");
 * \endcode
 *
 * Example with a string variable:
 *
 * \code
 * // -> retrieves global integer value "NumberOfCellsX" which is
 * // located in the parameter group "Grid"
 * std::string paramName = "Grid";
 * paramName += ".NumberOfCellsX";
 * GET_RUNTIME_PARAM_CSTRING(TypeTag, int, paramName.c_str());
 * \endcode
 */
#define GET_RUNTIME_PARAM_CSTRING(TypeTag, ParamType, ParamName) \
    ::Dumux::template getParam_UsingDeprecatedMacro<ParamType>(std::string(ParamName))

/*!
 * \ingroup Parameter
 * \brief Retrieve a runtime parameter which _does not_ have a default value taken from
 *        the Dumux property system.
 *
 * The third argument is group name.
 *
 * Example:
 *
 * \code
 * // -> retrieves global integer value "NumberOfCellsX" which is
 * // located in the parameter group "Grid"
 * GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NumberOfCellsX);
 * \endcode
 */
#define GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, ParamType, GroupName, ParamName) \
    ::Dumux::template getParam_UsingDeprecatedMacro<ParamType>(std::string(#GroupName) + "." + std::string(#ParamName))

/*!
 * \ingroup Parameter
 * \brief Retrieve a runtime parameter which _does not_ have a default value taken from
 *        the Dumux property system.
 *
 * The third argument is group name, which has to be a c-string.
 * This allows to use string variables as group name. The functionality of having variables as
 * group name is no problem when directly using the Dune::ParameterTree. For consistency with the
 * macro way of reading in the parameters this macro is necessary e.g. in the gridcreator to reach
 * a satisfying level of generality.
 *
 * Example with a temporary c-string:
 *
 * \code
 * // -> retrieves global integer value "NumberOfCellsX" which is
 * // located in the parameter group "Grid"
 * GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, int, "Grid", NumberOfCellsX);
 * \endcode
 *
 * Example with a string variable:
 *
 * \code
 * // -> retrieves global integer value "NumberOfCellsX" which is
 * // located in the parameter group "Grid"
 * std::string groupName = "Grid";
 * GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, int, groupName.c_str(), NumberOfCellsX);
 * \endcode
 */
#define GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, ParamType, GroupName, ParamName) \
    ::Dumux::template getParam_UsingDeprecatedMacro<ParamType>(std::string(GroupName) + "." + std::string(#ParamName))

namespace Dumux
{

//! The runtime parameter managing class
class Parameters {

public:

    //! Initialize the parameter tree singletons
    static void init(int argc, char **argv,
                     std::string parameterFileName = "",
                     void (*usage)(const char *, const std::string &) = [](const char *, const std::string &){})
    {
        const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

        // check whether the user wanted to see the help message
        for (int i = 1; i < argc; ++i)
        {
            if (std::string("--help") == argv[i] || std::string("-h") == argv[i])
            {
                // return usage message and return;
                if (mpiHelper.rank() == 0)
                    usage(argv[0], defaultUsageMessage(argv[0]));

                exit(0);
            }
        }

        // apply the default parameters
        defaultParameters(defaultParamTree());

        // parse paramters from the command line
        for (int i = 1; i < argc; ++i)
        {
            if (argv[i][0] != '-' && i == 1)
            {
                // try to pass first argument as parameter file
                parameterFileName = argv[1];
                break;
            }

            if (argv[i][0] != '-')
                DUNE_THROW(ParameterException, "-> Command line argument " << i << " (='" << argv[i] << "') is invalid. <-");

            if (i+1 == argc)
                DUNE_THROW(ParameterException, "-> No argument given for parameter '" << argv[i] << "'! <-");

            // check for the ParameterFile argument
            if (argv[i]+1 == std::string("ParameterFile")) // +1 removes the '-'
            {
                parameterFileName = argv[i+1];
                ++i;
            }

            // add all other options as key value pairs
            else
            {
                // read a -MyOpt VALUE option
                std::string paramName = argv[i]+1; // +1 removes the '-'
                std::string paramValue = argv[i+1];
                ++i; // In the case of '-MyOpt VALUE' each pair counts as two arguments

                // Put the key=value pair into the parameter tree
                paramTree()[paramName] = paramValue;
            }
        }

        // otherwise use the default name (executable name + .input)
        if (parameterFileName == "")
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
            // read parameters from the file without overwriting the command line params
            // because the command line arguments have precedence
            Dune::ParameterTreeParser::readINITree(parameterFileName,
                                                   paramTree(),
                                                   /*overwrite=*/false);
        }
        parameterFile.close();
    }

    //! returns the logging parameter tree recording which parameters are used during the simulation
    static const LoggingParameterTree& getTree()
    {
        static LoggingParameterTree tree(paramTree(), defaultParamTree());
        return tree;
    }

private:

    //! the actual internal parameter tree storing all user-specfied runtime parameters
    static Dune::ParameterTree& paramTree()
    {
        static Dune::ParameterTree tree;
        return tree;
    }

    //! the parameter tree storing the Dumux global defaults for some parameters
    static Dune::ParameterTree& defaultParamTree()
    {
        static Dune::ParameterTree tree;
        return tree;
    }

    //! This method puts all default arguments into the parameter tree
    //! we do this once per simulation on call to Parameters::init();
    static void defaultParameters(Dune::ParameterTree& params)
    {
        // parameters in the problem group
        params["Problem.EnableGravity"] = "true";

        // parameters in the newton group
        params["Newton.TargetSteps"] = "16";
    }
};

enum class ParamLookup {
    simple, tree
};

// a free function to get a parameter from the parameter tree singleton
// e.g. auto endTime = getParam<double>("TimeManager.TEnd");
template<typename T, typename... Args>
T getParam(Args&&... args)
{
    const auto& p = Parameters::getTree();
    return p.template get<T>(std::forward<Args>(args)... );
}

template<typename T, typename... Args>
DUNE_DEPRECATED_MSG("Using preprocessor MACROS for getting parameters is deprecated on next. Please use the new getParam method.")
T getParam_UsingDeprecatedMacro(Args&&... args)
{
    const auto& p = Parameters::getTree();
    return p.template get<T>(std::forward<Args>(args)... );
}

// a free function to report all parameters stating currently unused ones
void reportParams()
{ Parameters::getTree().reportAll(); }

} // namespace Dumux

#endif
