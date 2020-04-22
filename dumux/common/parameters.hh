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
#include <functional>

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
                     const Usage& usage)
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
     *         1) global defaults (see member function applyGlobalDefaults_)
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
        defaultParams(defaultParamTree_());
        applyGlobalDefaults_(defaultParamTree_());

        // parse paramters from the command line
        const auto commandLineArgs = parseCommandLine(argc, argv);
        mergeTree_(paramTree_(), commandLineArgs);

        // overwrite parameter file if one was specified on the command line
        parameterFileName = commandLineArgs.get<std::string>("ParameterFile", parameterFileName);

        // otherwise use the default name (executable name + .input)
        if (parameterFileName.empty())
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
            if (parameterFileName.empty())
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
                                               paramTree_(),
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
        params(paramTree_());
        // apply the default parameters
        defaultParams(defaultParamTree_());
        applyGlobalDefaults_(defaultParamTree_());
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
        params(paramTree_());

        // read parameters from the input file
        Dune::ParameterTreeParser::readINITree(parameterFileName, paramTree_(), inputFileOverwritesParams);

        // apply the default parameters
        defaultParams(defaultParamTree_());
        applyGlobalDefaults_(defaultParamTree_());
    }

    //! prints all used and unused parameters
    static void print()
    {
        getTree().reportAll();
    }

    //! Parse command line arguments into a parameter tree
    static Dune::ParameterTree parseCommandLine(int argc, char **argv)
    {
        Dune::ParameterTree commandLineArgs;
        for (int i = 1; i < argc; ++i)
        {
            if (argv[i][0] != '-' && i == 1)
            {
                // try to pass first argument as parameter file
                commandLineArgs["ParameterFile"] = argv[1];
                continue;
            }

            if (argv[i][0] != '-')
                DUNE_THROW(ParameterException, "-> Command line argument " << i << " (='" << argv[i] << "') is invalid. <-");

            if (i+1 == argc)
                DUNE_THROW(ParameterException, "-> No argument given for parameter '" << argv[i] << "'! <-");

            // check for the ParameterFile argument
            if (argv[i]+1 == std::string("ParameterFile")) // +1 removes the '-'
            {
                commandLineArgs["ParameterFile"] = argv[i+1];
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
                commandLineArgs[paramName] = paramValue;
            }
        }
        return commandLineArgs;
    }

    /*!
     * \brief Get the parameter tree
     *
     * The logging parameter tree recording which parameters are used during the simulation
     */
    static const LoggingParameterTree& getTree()
    {
        static LoggingParameterTree tree(paramTree_(), defaultParamTree_());
        return tree;
    }

private:
    //! the actual internal parameter tree storing all user-specfied runtime parameters
    static Dune::ParameterTree& paramTree_()
    {
        static Dune::ParameterTree tree;
        return tree;
    }

    //! the parameter tree storing the Dumux global defaults for some parameters
    static Dune::ParameterTree& defaultParamTree_()
    {
        static Dune::ParameterTree tree;
        return tree;
    }

    //! This method puts all default arguments into the parameter tree
    //! we do this once per simulation on call to Parameters::init();
    static void applyGlobalDefaults_(Dune::ParameterTree& params)
    {
        // global defaults
        Dune::ParameterTree defaultParams;

        // parameters in the implicit group
        defaultParams["Flux.UpwindWeight"] = "1.0";
        defaultParams["Implicit.EnableJacobianRecycling"] = "false";

        // parameters in the assembly group
        defaultParams["Assembly.NumericDifferenceMethod"] = "1";

        // parameters in the problem group
        defaultParams["Problem.EnableGravity"] = "true";
        defaultParams["Problem.EnableInertiaTerms"] = "true";

        // parameters in the Newton group
        // MinSteps = 2 makes Newton more robust if converge criterion is not perfect
        defaultParams["Newton.MinSteps"] = "2";
        defaultParams["Newton.MaxSteps"] = "18";
        defaultParams["Newton.TargetSteps"] = "10";
        defaultParams["Newton.UseLineSearch"] = "false";
        defaultParams["Newton.EnableChop"] = "false";
        defaultParams["Newton.EnableShiftCriterion"] = "true";
        defaultParams["Newton.MaxRelativeShift"] = "1e-8";
        defaultParams["Newton.EnableResidualCriterion"] = "false";
        defaultParams["Newton.ResidualReduction"] = "1e-5";
        defaultParams["Newton.EnableAbsoluteResidualCriterion"] = "false";
        defaultParams["Newton.MaxAbsoluteResidual"] = "1e-5";
        defaultParams["Newton.SatisfyResidualAndShiftCriterion"] = "false";
        defaultParams["Newton.EnablePartialReassembly"] = "false";

        // parameters in the time loop group
        defaultParams["TimeLoop.MaxTimeStepSize"] = "1e300";
        defaultParams["TimeLoop.MaxTimeStepDivisions"] = "10";

        // parameters in the vtk group
        defaultParams["Vtk.AddVelocity"] = "false";
        defaultParams["Vtk.AddProcessRank"] = "true";

        // parameters in the mpfa group
        defaultParams["MPFA.Q"] = "0.0";

        // merge the global default tree but do not overwrite if the parameter already exists
        mergeTree_(params, defaultParams, false);
    }

    //! merge source into target tree
    static void mergeTree_(Dune::ParameterTree& target, const Dune::ParameterTree& source, bool overwrite = true)
    { mergeTreeImpl_(target, source, overwrite, ""); }

    //! recursively merge all elements
    static void mergeTreeImpl_(Dune::ParameterTree& target, const Dune::ParameterTree& source, bool overwrite, const std::string& group)
    {
        const auto prefix = group.empty() ? "" : group + ".";
        for (const auto& key : source.getValueKeys())
            if (overwrite || !target.hasKey(key))
                target[prefix + key] = source[key];

        for (const auto& subKey : source.getSubKeys())
            mergeTreeImpl_(target, source.sub(subKey), overwrite, prefix + subKey);
    }
};

/*!
 * \ingroup Common
 * \brief A free function to get a parameter from the parameter tree singleton
 * \note \code auto endTime = getParam<double>("TimeManager.TEnd"); \endcode
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
template<typename T, typename... Args>
T getParam(Args&&... args)
{ return Parameters::getTree().template get<T>(std::forward<Args>(args)... ); }

/*!
 * \ingroup Common
 * \brief A free function to get a parameter from the parameter tree singleton with a model group
 * \note \code  auto endTime = getParamFromGroup<double>("FreeFlow", "TimeManager.TEnd"); \endcode
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
template<typename T, typename... Args>
T getParamFromGroup(Args&&... args)
{ return Parameters::getTree().template getFromGroup<T>(std::forward<Args>(args)... ); }

/*!
 * \ingroup Common
 * \brief Check whether a key exists in the parameter tree
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
inline bool hasParam(const std::string& param)
{ return Parameters::getTree().hasKey(param); }

/*!
 * \ingroup Common
 * \brief Check whether a key exists in the parameter tree with a model group prefix
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
inline bool hasParamInGroup(const std::string& paramGroup, const std::string& param)
{ return Parameters::getTree().hasKeyInGroup(param, paramGroup); }

/*!
 * \ingroup Common
 * \brief Get a list of sub groups from the parameter tree sorted by relevance
 * \return A vector of fully qualified subGroup names sorted by descending relevance.
 */
inline std::vector<std::string> getParamSubGroups(const std::string& subGroupName, const std::string& paramGroup)
{ return Parameters::getTree().getSubGroups(subGroupName, paramGroup); }

} // end namespace Dumux

#endif
