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
 * \brief A parameter tree that logs which parameters have been used
 */
#ifndef DUMUX_LOGGING_PARAMETER_TREE_HH
#define DUMUX_LOGGING_PARAMETER_TREE_HH

#include <iomanip>
#include <iostream>
#include <string>

#include <dune/common/parametertree.hh>

namespace Dumux
{

/*!
 * \ingroup Common
 * \brief A parameter tree that logs which parameters have been used
 */
class LoggingParameterTree
{

public:
    /*
     * \brief A logging parameter tree is always attached to an existingparameter tree
     */
    LoggingParameterTree() = delete;

    /*
     * \brief Create LoggingParameterTree from ParameterTree
     */
    LoggingParameterTree(const Dune::ParameterTree& params, const Dune::ParameterTree& defaultParams)
    : params_(params), defaultParams_(defaultParams) {}

    /** \brief test for key
     *
     * Tests whether given key exists.
     *
     * \param key key name
     * \return true if key exists in structure, otherwise false
     */
    bool hasKey(const std::string& key) const
    { return params_.hasKey(key); }


    /** \brief print the hierarchical parameter tree to stream
     *
     * \param stream the output stream to print to
     */
    void report(std::ostream& stream = std::cout) const
    { params_.report(stream); }

    /** \brief print distinct substructure to stream
     *
     * Prints all entries with given prefix.
     *
     * \param stream Stream to print to
     * \param prefix for key and substructure names
     */
    void reportAll(std::ostream& stream = std::cout) const
    {
        stream << "\n# Runtime-specified parameters used:" << std::endl;
        usedRuntimeParams_.report(stream);

        const auto unusedParams = getUnusedKeys();
        if (!unusedParams.empty())
        {
            stream << "\n# Unused parameters:" << std::endl;
            for (const auto& key : unusedParams)
                stream << key << " = \"" << params_[key] << "\"" << std::endl;
        }
    }


    /** \brief get value as string
     *
     * Returns pure string value for given key.
     *
     * \param key key name
     * \param defaultValue default if key does not exist
     * \return value as string
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    std::string get(const std::string& key, const std::string& defaultValue) const
    {
        if (params_.hasKey(key))
        {
            // log that we used this parameter
            usedRuntimeParams_[key] = params_[key];
            return params_[key];
        }

        return defaultValue;
    }

    /** \brief get value as string
     *
     * Returns pure string value for given key.
     *
     * \todo This is a hack so get("my_key", "xyz") compiles
     * (without this method "xyz" resolves to bool instead of std::string)
     * \param key key name
     * \param defaultValue default if key does not exist
     * \return value as string
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    std::string get(const std::string& key, const char* defaultValue) const
    {
        if (params_.hasKey(key))
        {
            // log that we used this parameter
            usedRuntimeParams_[key] = params_[key];
            return params_[key];
        }

        return defaultValue;
    }


    /** \brief get value converted to a certain type
     *
     * Returns value as type T for given key.
     *
     * \tparam T type of returned value.
     * \param key key name
     * \param defaultValue default if key does not exist
     * \return value converted to T
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    template<typename T>
    T get(const std::string& key, const T& defaultValue) const
    {
        if (params_.hasKey(key))
        {
            // log that we used this parameter
            usedRuntimeParams_[key] = params_[key];
            return params_.template get<T>(key);
        }

        return defaultValue;
    }

    /** \brief Get value
     *
     * \tparam T Type of the value
     * \param key Key name
     * \throws RangeError if key does not exist
     * \throws NotImplemented Type is not supported
     * \return value as T
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    template <class T>
    T get(const std::string& key) const
    {
        if (params_.hasKey(key))
        {
            // log that we used this parameter
            usedRuntimeParams_[key] = params_[key];
        }

        else if(defaultParams_.hasKey(key))
        {
            // use the default
            return defaultParams_.template get<T>(key);
        }

        return params_.template get<T>(key);
    }

    /** \brief Find the keys that haven't been used yet
     *
     * \retuns unusedParams Container storing unused keys
     * \note Useful for debugging purposes
     */
    std::vector<std::string> getUnusedKeys() const
    {
        std::vector<std::string> unusedParams;
        findUnusedKeys(params_, unusedParams);
        return unusedParams;
    }

private:
    /** \brief Find the keys that haven't been used yet recursively
     *
     * \param tree The tree to look in for unused keys
     * \param unusedParams Container to store unused keys
     * \param prefix the prefix attached to the key
     */
    void findUnusedKeys(const Dune::ParameterTree& tree,
                        std::vector<std::string>& unusedParams,
                        const std::string& prefix = "") const
    {
        // loop over all keys of the current tree
        // store keys which were not accessed
        const auto& keys = tree.getValueKeys();
        for (const auto& key : keys)
            if (!usedRuntimeParams_.hasKey(prefix + key))
                unusedParams.push_back(prefix + key);

        // recursively loop over all subtrees
        const auto& subTreeKeys = tree.getSubKeys();
        for (const auto& key : subTreeKeys)
            findUnusedKeys(tree.sub(key), unusedParams, prefix + key + ".");
    }

    const Dune::ParameterTree& params_;
    const Dune::ParameterTree& defaultParams_;
    mutable Dune::ParameterTree usedRuntimeParams_;
};

} // end namespace Dumux

#endif
