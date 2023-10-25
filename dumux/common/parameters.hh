// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parameter
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

#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/loggingparametertree.hh>

namespace Dumux {

/*!
 * \ingroup Parameter
 * \brief Parameter class managing runtime input parameters
 * \todo Doc me!
 */
class Parameters {

    using DefaultParams = std::function<void (Dune::ParameterTree&)>;
    using Usage = std::function<void (const char *, const std::string &)>;

public:

    //! Initialize the parameter tree singletons
    static void init(int argc, char **argv, const Usage& usage);

    //! Initialize the parameter tree singletons
    static void init(int argc, char **argv,
                    std::string parameterFileName,
                    const Usage& usage = [](const char *, const std::string &){});

    //! Initialize the parameter tree singletons
    static void init(int argc, char **argv,
                     const DefaultParams& defaultParams,
                     const Usage& usage);

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
                     const Usage& usage = [](const char *, const std::string &){});

    /*!
     * \brief Initialize the parameter tree
     * \param params a function that sets parameters of the runtime parameter tree
     * \param defaultParams a function that sets parameters of the default runtim parameter tree
     * \note if a parameter is looked up without explicitly providing a default, the
     *       default tree is consulted if the parameter could not be found in the parameter tree
     */
    static void init(const DefaultParams&  params = [] (Dune::ParameterTree&) {},
                     const DefaultParams& defaultParams = [] (Dune::ParameterTree&) {});

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
                     const DefaultParams& defaultParams = [] (Dune::ParameterTree&) {});

    //! prints all used and unused parameters
    static void print();

    //! Parse command line arguments into a parameter tree
    static Dune::ParameterTree parseCommandLine(int argc, char **argv);

    /*!
     * \brief Get the parameter tree
     *
     * The logging parameter tree recording which parameters are used during the simulation
     */
    static const LoggingParameterTree& getTree();

private:
    //! the actual internal parameter tree storing all user-specfied runtime parameters
    static Dune::ParameterTree& paramTree_();

    //! the parameter tree storing the Dumux global defaults for some parameters
    static Dune::ParameterTree& defaultParamTree_();

    //! This method puts all default arguments into the parameter tree
    //! we do this once per simulation on call to Parameters::init();
    static void applyGlobalDefaults_(Dune::ParameterTree& params);

    //! merge source into target tree
    static void mergeTree_(Dune::ParameterTree& target, const Dune::ParameterTree& source, bool overwrite = true);

    //! recursively merge all elements
    static void mergeTreeImpl_(Dune::ParameterTree& target, const Dune::ParameterTree& source, bool overwrite, const std::string& group);
};

/*!
 * \ingroup Parameter
 * \brief A free function to get a parameter from the parameter tree singleton
 * \note \code auto endTime = getParam<double>("TimeManager.TEnd"); \endcode
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
template<typename T = std::string, typename... Args>
T getParam(Args&&... args)
{ return Parameters::getTree().template get<T>(std::forward<Args>(args)... ); }

/*!
 * \ingroup Parameter
 * \brief A free function to get a parameter from the parameter tree singleton with a model group
 * \note \code  auto endTime = getParamFromGroup<double>("FreeFlow", "TimeManager.TEnd"); \endcode
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
template<typename T = std::string, typename... Args>
T getParamFromGroup(Args&&... args)
{ return Parameters::getTree().template getFromGroup<T>(std::forward<Args>(args)... ); }

/*!
 * \ingroup Parameter
 * \brief Check whether a key exists in the parameter tree
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
inline bool hasParam(const std::string& param)
{ return Parameters::getTree().hasKey(param); }

/*!
 * \ingroup Parameter
 * \brief Check whether a key exists in the parameter tree with a model group prefix
 * \note Once this has been called the first time, you cannot modify the parameter tree anymore
 */
inline bool hasParamInGroup(const std::string& paramGroup, const std::string& param)
{ return Parameters::getTree().hasKeyInGroup(param, paramGroup); }

/*!
 * \ingroup Parameter
 * \brief Get a list of sub groups from the parameter tree sorted by relevance
 * \return A vector of fully qualified subGroup names sorted by descending relevance.
 */
inline std::vector<std::string> getParamSubGroups(const std::string& subGroupName, const std::string& paramGroup)
{ return Parameters::getTree().getSubGroups(subGroupName, paramGroup); }

} // end namespace Dumux

#endif
