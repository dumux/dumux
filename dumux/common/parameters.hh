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
 * \brief The infrastructure to retrieve run-time parameters from Dune::ParameterTrees.
 */
#ifndef DUMUX_PARAMETERS_HH
#define DUMUX_PARAMETERS_HH

#include <iostream>
#include <list>
#include <sstream>
#include <unordered_map>
#include <fstream>

#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/loggingparametertree.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Parameter class managing runtime input parameters
 * \todo Doc me!
 */
class Parameters {

    using DefaultParams = std::function<void (Dune::ParameterTree&)>;
    using Usage = std::function<void (const char *, const std::string &)>;

public:

    //! Initialize the parameter tree singletons
    static void init(int argc, char **argv, const Usage& usage)
    {
        init(argc, argv, [] (Dune::ParameterTree&) {}, "", usage);
    }

    //! Initialize the parameter tree singletons
    static void init(int argc, char **argv,
                    std::string parameterFileName,
                    const Usage& usage = [](const char *, const std::string &){})
    {
        init(argc, argv, [] (Dune::ParameterTree&) {}, parameterFileName, usage);
    }

    //! Initialize the parameter tree singletons
    static void init(int argc, char **argv,
                     const DefaultParams& defaultParams,
                     const Usage& usage = [](const char *, const std::string &){})
    {
        init(argc, argv, defaultParams, "", usage);
    }

    /*!
     * \brief Initialize the parameter tree
     * \param argc number of command line argument (forwarded from main)
     * \param argv command line argument (forwarded from main)
     * \param defaultParams a function that sets parameters of the default runtime parameter tree
     * \param parameterFileName the file name of the input file
     * \param usage the usage function to print if the help option was passed on the command line
     * \note the default parameter tree is initialized in the following way
     *         1) global defaults (see member function globalDefaultParameters)
     *         2) user provided defaults (overwrite global defaults)
     *       the parameter tree is initialized in the following way
     *         1) parameters from the input file
     *         2) parameters from the command line (overwrite input file parameters)
     * \note if a parameter is looked up without explicitly providing a default, the
     *       default tree is consulted if the parameter could not be found in the parameter tree
     */
    static void init(int argc, char **argv,
                     const DefaultParams& defaultParams = [] (Dune::ParameterTree&) {},
                     std::string parameterFileName = "",
                     const Usage& usage = [](const char *, const std::string &){})
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
        globalDefaultParameters(defaultParamTree());
        defaultParams(defaultParamTree());

        // parse paramters from the command line
        parameterFileName = parseCommandLineArguments(argc, argv, parameterFileName);

        // otherwise use the default name (executable name + .input)
        if (parameterFileName == "")
        {
            parameterFileName = [&](){
                std::string defaultName = std::string(argv[0]) + ".input";
                std::ifstream pFile(defaultName.c_str());
                if (pFile.is_open())
                    return defaultName;

                defaultName = "params.input";
                pFile = std::ifstream(defaultName.c_str());
                if (pFile.is_open())
                    return defaultName;
                else
                    return std::string("");
            }();

            // if no parameter file was given and also no default names where found, continue without
            if (parameterFileName == "")
            {
                if (mpiHelper.size() > 1)
                    std::cout << "Rank " << mpiHelper.rank() << ": ";
                std::cout << "No parameter file found. Continuing without parameter file.\n";

                return;
            }
            else
            {
                if (mpiHelper.size() > 1)
                    std::cout << "Rank " << mpiHelper.rank() << ": ";
                std::cout << "No parameter file given. "
                          << "Defaulting to '"
                          << parameterFileName
                          << "' for input file.\n";
            }
        }

        if (mpiHelper.size() > 1)
            std::cout << "Rank " << mpiHelper.rank() << ": ";
        std::cout << "Reading parameters from file " << parameterFileName << ".\n";

        // read parameters from the file without overwriting the command line params
        // because the command line arguments have precedence
        // let Dune do the error checking if the file exists
        Dune::ParameterTreeParser::readINITree(parameterFileName,
                                               paramTree(),
                                               /*overwrite=*/false);
    }

    /*!
     * \brief Initialize the parameter tree
     * \param params a function that sets parameters of the runtime parameter tree
     * \param defaultParams a function that sets parameters of the default runtim parameter tree
     * \note if a parameter is looked up without explicitly providing a default, the
     *       default tree is consulted if the parameter could not be found in the parameter tree
     */
    static void init(const DefaultParams&  params = [] (Dune::ParameterTree&) {},
                     const DefaultParams& defaultParams = [] (Dune::ParameterTree&) {})
    {
        // apply the parameters
        params(paramTree());
        // apply the default parameters
        globalDefaultParameters(defaultParamTree());
        defaultParams(defaultParamTree());
    }

    /*!
     * \brief Initialize the parameter tree
     * \param parameterFileName an input parameter file name
     * \param params a parameter tree with runtime parameters
     * \param inputFileOverwritesParams if set to true (default) the parameters from the input file have precedence,
     *                                  if set to false the input the parameters provided via params have precedence
     * \param defaultParams a parameter tree with default parameters
     * \note the params function overwrites
     * \note if a parameter is looked up without explicitly providing a default, the
     *       default tree is consulted if the parameter could not be found in the parameter tree
     */
    static void init(const std::string& parameterFileName,
                     const DefaultParams& params = [] (Dune::ParameterTree&) {},
                     bool inputFileOverwritesParams = true,
                     const DefaultParams& defaultParams = [] (Dune::ParameterTree&) {})
    {
        // apply the parameters
        params(paramTree());

        // read parameters from the input file
        Dune::ParameterTreeParser::readINITree(parameterFileName, paramTree(), inputFileOverwritesParams);

        // apply the default parameters
        globalDefaultParameters(defaultParamTree());
        defaultParams(defaultParamTree());
    }

    //! \brief parse the arguments given on the command line
    //! \returns the parameterFileName if one was given otherwise returns empty string
    static std::string parseCommandLineArguments(int argc, char **argv,
                                                 std::string parameterFileName = "")
    {
        for (int i = 1; i < argc; ++i)
        {
            if (argv[i][0] != '-' && i == 1)
            {
                // try to pass first argument as parameter file
                parameterFileName = argv[1];
                continue;
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
        return parameterFileName;
    }

    //! prints all used and unused parameters
    static void print()
    {
        getTree().reportAll();
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
    static void globalDefaultParameters(Dune::ParameterTree& params)
    {
        // parameters in the implicit group
        params["Flux.UpwindWeight"] = "1.0";
        params["Implicit.EnableJacobianRecycling"] = "false";

        // parameters in the assembly group
        params["Assembly.NumericDifferenceMethod"] = "1";

        // parameters in the linear solver group
        params["LinearSolver.GMResRestart"] = "10";
        params["LinearSolver.MaxIterations"] = "250";
        params["LinearSolver.PreconditionerIterations"] = "1";
        params["LinearSolver.PreconditionerRelaxation"] = "1.0";
        params["LinearSolver.ResidualReduction"] = "1e-13";
        params["LinearSolver.Verbosity"] = "0";

        // parameters in the problem group
        params["Problem.EnableGravity"] = "true";
        params["Problem.EnableInertiaTerms"] = "true";

        // parameters in the Newton group
        params["Newton.MinSteps"] = "0";
        params["Newton.MaxSteps"] = "18";
        params["Newton.TargetSteps"] = "10";
        params["Newton.UseLineSearch"] = "false";
        params["Newton.EnableChop"] = "false";
        params["Newton.EnableShiftCriterion"] = "true";
        params["Newton.MaxRelativeShift"] = "1e-8";
        params["Newton.EnableResidualCriterion"] = "false";
        params["Newton.ResidualReduction"] = "1e-5";
        params["Newton.EnableAbsoluteResidualCriterion"] = "false";
        params["Newton.MaxAbsoluteResidual"] = "1e-5";
        params["Newton.SatisfyResidualAndShiftCriterion"] = "false";
        params["Newton.EnablePartialReassembly"] = "false";

        // parameters in the time loop group
        params["TimeLoop.MaxTimeStepSize"] = "1e300";
        params["TimeLoop.MaxTimeStepDivisions"] = "10";

        // parameters in the vtk group
        params["Vtk.AddVelocity"] = "false";
        params["Vtk.AddProcessRank"] = "true";

        // parameters in the mpfa group
        params["Mpfa.Q"] = "0.0";
    }
};

/*!
 * \ingroup Common
 * \brief a free function to set model- or problem-specific default parameters
 */
void setParam(Dune::ParameterTree& params,
              const std::string& group,
              const std::string& key,
              const std::string& value)
{
    if(group == "")
        params[key] = value;
    else
        params[group + "." + key] = value;
}

/*!
 * \ingroup Common
 * \brief A free function to get a parameter from the parameter tree singleton
 * \note \code auto endTime = getParam<double>("TimeManager.TEnd"); \endcode
 */
template<typename T, typename... Args>
T getParam(Args&&... args)
{
    const auto& p = Parameters::getTree();
    return p.template get<T>(std::forward<Args>(args)... );
}

/*!
 * \ingroup Common
 * \brief A free function to get a parameter from the parameter tree singleton with a model group
 * \note \code  auto endTime = getParamFromGroup<double>("FreeFlow", "TimeManager.TEnd"); \endcode
 */
template<typename T, typename... Args>
T getParamFromGroup(Args&&... args)
{
    const auto& p = Parameters::getTree();
    return p.template getFromGroup<T>(std::forward<Args>(args)... );
}

/*!
 * \ingroup Common
 * \brief Check whether a key exists in the parameter tree
 */
bool hasParam(const std::string& param)
{
    const auto& p = Parameters::getTree();
    return p.hasKey(param);
}

/*!
 * \ingroup Common
 * \brief Check whether a key exists in the parameter tree with a model group prefix
 */
template<typename... Args>
bool hasParamInGroup(const std::string& paramGroup, const std::string& param)
{
    const auto& p = Parameters::getTree();
    if (paramGroup == "")
        return p.hasKey(param);
    else
        return p.hasKey(paramGroup + "." + param);
}

DUNE_DEPRECATED_MSG("haveParam is deprecated, please use hasParam instead.")
bool haveParam(const std::string& param)
{ return hasParam(param); }

template<typename... Args>
DUNE_DEPRECATED_MSG("haveParamInGroup is deprecated, please use hasParamInGroup instead.")
bool haveParamInGroup(const std::string& paramGroup, const std::string& param)
{ return hasParamInGroup(paramGroup, param); }

} // namespace Dumux

#endif
