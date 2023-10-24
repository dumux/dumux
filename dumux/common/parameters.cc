// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief The infrastructure to retrieve run-time parameters from Dune::ParameterTrees.
 */

#include <config.h>

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
#include <dumux/common/parameters.hh>
#include <dumux/common/loggingparametertree.hh>

namespace Dumux {

// Initialize the parameter tree singletons
void Parameters::init(int argc, char **argv, const Usage& usage)
{
    init(argc, argv, [] (Dune::ParameterTree&) {}, "", usage);
}

// Initialize the parameter tree singletons
void Parameters::init(int argc, char **argv,
                      std::string parameterFileName,
                      const Usage& usage)
{
    init(argc, argv, [] (Dune::ParameterTree&) {}, parameterFileName, usage);
}

// Initialize the parameter tree singletons
void Parameters::init(int argc, char **argv,
                      const DefaultParams& defaultParams,
                      const Usage& usage)
{
    init(argc, argv, defaultParams, "", usage);
}

// Initialize the parameter tree
void Parameters::init(int argc, char **argv,
                      const DefaultParams& defaultParams,
                      std::string parameterFileName,
                      const Usage& usage)
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

    // parse parameters from the command line
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

// Initialize the parameter tree
void Parameters::init(const DefaultParams& params,
                      const DefaultParams& defaultParams)
{
    // apply the parameters
    params(paramTree_());
    // apply the default parameters
    defaultParams(defaultParamTree_());
    applyGlobalDefaults_(defaultParamTree_());
}

// Initialize the parameter tree
void Parameters::init(const std::string& parameterFileName,
                      const DefaultParams& params,
                      bool inputFileOverwritesParams,
                      const DefaultParams& defaultParams)
{
    // apply the parameters
    params(paramTree_());

    // read parameters from the input file
    Dune::ParameterTreeParser::readINITree(parameterFileName, paramTree_(), inputFileOverwritesParams);

    // apply the default parameters
    defaultParams(defaultParamTree_());
    applyGlobalDefaults_(defaultParamTree_());
}

// prints all used and unused parameters
void Parameters::print()
{
    getTree().reportAll();
}

// Parse command line arguments into a parameter tree
Dune::ParameterTree Parameters::parseCommandLine(int argc, char **argv)
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

// get the parameter tree singleton
const LoggingParameterTree& Parameters::getTree()
{
    static LoggingParameterTree tree(paramTree_(), defaultParamTree_());
    return tree;
}

// the actual internal parameter tree storing all user-specfied runtime parameters
Dune::ParameterTree& Parameters::paramTree_()
{
    static Dune::ParameterTree tree;
    return tree;
}

// the parameter tree storing the Dumux global defaults for some parameters
Dune::ParameterTree& Parameters::defaultParamTree_()
{
    static Dune::ParameterTree tree;
    return tree;
}

void Parameters::applyGlobalDefaults_(Dune::ParameterTree& params)
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

// merge source into target tree
void Parameters::mergeTree_(Dune::ParameterTree& target, const Dune::ParameterTree& source, bool overwrite)
{ mergeTreeImpl_(target, source, overwrite, ""); }

// recursively merge all elements
void Parameters::mergeTreeImpl_(Dune::ParameterTree& target, const Dune::ParameterTree& source, bool overwrite, const std::string& group)
{
    const auto prefix = group.empty() ? "" : group + ".";
    for (const auto& key : source.getValueKeys())
        if (overwrite || !target.hasKey(key))
            target[prefix + key] = source[key];

    for (const auto& subKey : source.getSubKeys())
        mergeTreeImpl_(target, source.sub(subKey), overwrite, prefix + subKey);
}

} // end namespace Dumux
