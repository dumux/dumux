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
 * \brief A parameter tree that logs which parameters have been used
 */
#ifndef DUMUX_LOGGING_PARAMETER_TREE_HH
#define DUMUX_LOGGING_PARAMETER_TREE_HH

#include <iomanip>
#include <iostream>
#include <string>

#include <dune/common/parametertree.hh>
#include <dumux/common/exceptions.hh>

namespace Dumux {

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
     * \note This ignores defaults. Hence, if the
     *       the specified key only exists in the defaults, this
     *       function returns false
     *
     * \param key key name
     * \return true if key exists in structure, otherwise false
     */
    bool hasKey(const std::string& key) const
    { return params_.hasKey(key); }

    /** \brief test for key in group
     *
     * Tests whether given key exists in a group.
     * Given a group this function starts to look from the back
     *       for dots. In G1.G2.G3 the function first looks if the key
     *       "G3.Key" exists, then "G2.Key", ...
     *
     * \note This ignores defaults. Hence, if the
     *       the specified key only exists in the defaults, this
     *       function returns false
     *
     * \param key key name
     * \param groupPrefix the group prefix name
     * \return true if key exists in structure, otherwise false
     */
    bool hasKeyInGroup(const std::string& key,
                       const std::string& groupPrefix) const
    {
        if (groupPrefix.empty())
            return hasKey(key);

        if (hasKey(key))
            return true;

        auto compoundKey = groupPrefix + "." + key;
        if (params_.hasKey(compoundKey))
            return true;

        compoundKey = findKeyInGroup(params_, key, groupPrefix);
        if (compoundKey != "")
            return true;

        return false;
    }

    /** \brief obtain a vector of all full group names for a specified subgroup name
     *
     * Example:
     * ------------
     * For the parameter tree
     *
     * [G1]
     * MyParam1 = 1
     * [G2.G1]
     * MyParam2 = 2
     * [G3.G2.G1]
     * MyParam3 = 3
     *
     * and groupPrefix="G3.G2" and subGroupName="G1"
     * this returns a vector with the entries {"G3.G2.G1", "G2.G1", "G1"}.
     * If groupPrefix = "G2", it returns {"G2.G1", "G1"}.
     * If groupPrefix = "" the returned vector has size 1 (containing subGroupName),
     * or size 0 if the subgroup does not exist in the parameter tree.
     *
     * \param subGroupName the sub group to look for
     * \param groupPrefix the group prefix name (potentially prefixing the subgroup)
     * \return a vector of fully qualified groups ordered by decreasing tree depth
     */
    std::vector<std::string> getSubGroups(const std::string& subGroupName,
                                          std::string groupPrefix) const
    {
        std::vector<std::string> groupNames;

        if (!groupPrefix.empty())
        {
            auto compoundGroup = groupPrefix + "." + subGroupName;
            for (std::string::size_type dotPos = 0; dotPos != std::string::npos; dotPos = groupPrefix.rfind("."))
            {
                if (params_.hasSub(compoundGroup) || defaultParams_.hasSub(compoundGroup))
                    groupNames.push_back(compoundGroup);

                groupPrefix = groupPrefix.substr(0, dotPos);
                compoundGroup = groupPrefix + "." + subGroupName;
            }
        }

        if (params_.hasSub(subGroupName) || defaultParams_.hasSub(subGroupName))
            groupNames.push_back(subGroupName);

        return groupNames;
    }

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
     */
    void reportAll(std::ostream& stream = std::cout) const
    {
        stream << "\n# Runtime-specified parameters used:" << std::endl;
        usedRuntimeParams_.report(stream);

        stream << "\n# Global default parameters used:" << std::endl;
        usedDefaultParams_.report(stream);

        const auto unusedParams = getUnusedKeys();
        if (!unusedParams.empty())
        {
            stream << "\n# Unused parameters:" << std::endl;
            for (const auto& key : unusedParams)
                stream << key << " = \"" << params_[key] << "\"" << std::endl;
        }
    }

    /*! \brief Do a backwards hierarchical search for a key in a group
     *
     * \note Given a group this function starts to look from the back
     *       for dots. In G1.G2.G3 the function first looks if the key
     *       "G3.Key" exists, then "G2.Key", ...
     *       The first compound key that is found is returned. If no
     *       compound key is found an empty string is returned.
     * \param tree The tree to look in for keys
     * \param key The key
     * \param groupPrefix the group prefix attached to the key
     */
    std::string findKeyInGroup(const Dune::ParameterTree& tree,
                               const std::string& key,
                               const std::string& groupPrefix) const
    {
        // search backwards until key is found
        std::string prefix = groupPrefix;
        auto dot = prefix.rfind(".");
        while (dot != std::string::npos)
        {
            prefix = prefix.substr(0, dot);
            std::string compoundKey = prefix + "." + key;

            if (tree.hasKey(compoundKey))
                return compoundKey;

            // look for the next dot in the current prefix
            dot = prefix.rfind(".");

        }

        // the key was not found
        return "";
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
            const auto returnValue = params_[key];
            usedRuntimeParams_[key] = returnValue;
            return returnValue;
        }

        return defaultValue;
    }

    /** \brief get value as string, preferably from the sub-tree corresponding
     *         to a given prefix. The sub-tree is searched backwards for the parameter
     *         until its "first" occurrence.
     *
     * Returns pure string value for given key.
     *
     * \param groupPrefix The prefix of the sub tree the search should start in
     * \param key key name
     * \param defaultValue default if key does not exist
     * \return value as string
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    std::string getFromGroup(const std::string& groupPrefix,
                             const std::string& key,
                             const std::string& defaultValue) const
    {
        if (groupPrefix.empty())
            return get(key, defaultValue);

        // first, look for the compound key
        std::string compoundKey = groupPrefix + "." + key;
        if (params_.hasKey(compoundKey))
        {
            // log that we used this parameter
            const auto returnValue = params_[compoundKey];
            usedRuntimeParams_[compoundKey] = returnValue;
            return returnValue;
        }

        // search backwards until key is found
        compoundKey = findKeyInGroup(params_, key, groupPrefix);
        if (compoundKey != "")
        {
            // log that we used this parameter
            const auto returnValue = params_[compoundKey];
            usedRuntimeParams_[compoundKey] = returnValue;
            return returnValue;
        }

        // finally, look for the key without prefix
        return get(key, defaultValue);
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
        const std::string dv = defaultValue;
        return get(key, dv);
    }

    /** \brief get value as string, preferably from the sub-tree corresponding
     *         to a given prefix. The sub-tree is searched for the parameter
     *         recursively until its "first" occurrence.
     *
     * Returns pure string value for given key.
     *
     * \param groupPrefix The prefix of the sub tree the search should start in
     * \param key key name
     * \param defaultValue default if key does not exist
     * \return value as string
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    std::string getFromGroup(const std::string& groupPrefix,
                             const std::string& key,
                             const char* defaultValue) const
    {
        const std::string dv = defaultValue;
        return getFromGroup(groupPrefix, key, dv);
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

    /** \brief get value as string, preferably from the sub-tree corresponding
     *         to a given prefix. The sub-tree is searched for the parameter
     *         recursively until its "first" occurrence.
     *
     * Returns pure string value for given key.
     *
     * \param groupPrefix The prefix of the sub tree the search should start in
     * \param key key name
     * \param defaultValue default if key does not exist
     * \return value as string
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    template<typename T>
    T getFromGroup(const std::string& groupPrefix,
                   const std::string& key,
                   const T& defaultValue) const
    {
        if (groupPrefix.empty())
            return get<T>(key, defaultValue);

        // first, look for the compound key
        std::string compoundKey = groupPrefix + "." + key;
        if (params_.hasKey(compoundKey))
        {
            // log that we used this parameter
            usedRuntimeParams_[compoundKey] = params_[compoundKey];
            return params_.template get<T>(compoundKey);
        }

        // search backwards until key is found
        compoundKey = findKeyInGroup(params_, key, groupPrefix);
        if (compoundKey != "")
        {
            // log that we used this parameter
            usedRuntimeParams_[compoundKey] = params_[compoundKey];
            return params_.template get<T>(compoundKey);
        }

        // finally, look for the key without prefix
        return get<T>(key, defaultValue);
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
            return params_.template get<T>(key);
        }

        else if(defaultParams_.hasKey(key))
        {
            // use the default
            usedDefaultParams_[key] = defaultParams_[key];
            return defaultParams_.template get<T>(key);
        }

        DUNE_THROW(Dumux::ParameterException, "Key " << key << " not found in the parameter tree");
    }

    /** \brief get value as string, preferably from the sub-tree corresponding
     *         to a given prefix. The sub-tree is searched for the parameter
     *         recursively until its "first" occurrence.
     *
     * Returns pure string value for given key.
     *
     * \param groupPrefix The prefix of the sub tree the search should start in
     * \param key key name
     * \return value as string
     * \note This might be quite slow so call it only once,
     *       e.g. on initialization of the class configured by runtime parameters
     */
    template<typename T>
    T getFromGroup(const std::string& groupPrefix,
                   const std::string& key) const
    {
        if (groupPrefix.empty())
            return get<T>(key);

        // first, look for the compound key
        std::string compoundKey = groupPrefix + "." + key;
        if (params_.hasKey(compoundKey))
        {
            // log that we used this parameter
            usedRuntimeParams_[compoundKey] = params_[compoundKey];
            return params_.template get<T>(compoundKey);
        }

        // search backwards until key is found
        compoundKey = findKeyInGroup(params_, key, groupPrefix);
        if (compoundKey != "")
        {
            // log that we used this parameter
            usedRuntimeParams_[compoundKey] = params_[compoundKey];
            return params_.template get<T>(compoundKey);
        }

        // reset the compoundKey
        compoundKey = groupPrefix + "." + key;

        // if the backward search did not succeed, try the bare key without any prefix
        if (params_.hasKey(key))
        {
            // log that we used this parameter
            usedRuntimeParams_[key] = params_[key];
            return params_.template get<T>(key);
        }

        // if this did not work, repeat the procedure using the default parameters
        else if(defaultParams_.hasKey(compoundKey))
        {
            // use the default
            usedDefaultParams_[compoundKey] = defaultParams_[compoundKey];
            return defaultParams_.template get<T>(compoundKey);
        }

        else
        {
            // search backwards until key is found
            compoundKey = findKeyInGroup(defaultParams_, key, groupPrefix);
            if (compoundKey != "")
            {
                // log that we used this parameter
                usedDefaultParams_[compoundKey] = defaultParams_[compoundKey];
                return defaultParams_.template get<T>(compoundKey);
            }

            if(defaultParams_.hasKey(key))
            {
                // use the default
                usedDefaultParams_[key] = defaultParams_[key];
                return defaultParams_.template get<T>(key);
            }

            DUNE_THROW(Dumux::ParameterException, "Key " << key << " not found in the parameter tree");
        }
    }

    /** \brief Find the keys that haven't been used yet
     *
     * \return unusedParams Container storing unused keys
     * \note Useful for debugging purposes
     */
    std::vector<std::string> getUnusedKeys() const
    {
        std::vector<std::string> unusedParams;
        findUnusedKeys_(params_, unusedParams);
        return unusedParams;
    }

private:
    /** \brief Find the keys that haven't been used yet recursively
     *
     * \param tree The tree to look in for unused keys
     * \param unusedParams Container to store unused keys
     * \param prefix the prefix attached to the key
     */
    void findUnusedKeys_(const Dune::ParameterTree& tree,
                         std::vector<std::string>& unusedParams,
                         const std::string& prefix = "") const
    {
        // loop over all keys of the current tree
        // store keys which were not accessed
        const auto& keys = tree.getValueKeys();
        for (const auto& key : keys)
            if (key != "ParameterFile" && !usedRuntimeParams_.hasKey(prefix + key))
                unusedParams.push_back(prefix + key);

        // recursively loop over all subtrees
        const auto& subTreeKeys = tree.getSubKeys();
        for (const auto& key : subTreeKeys)
            findUnusedKeys_(tree.sub(key), unusedParams, prefix + key + ".");
    }

    const Dune::ParameterTree& params_;
    const Dune::ParameterTree& defaultParams_;

    // logging caches
    mutable Dune::ParameterTree usedRuntimeParams_;
    mutable Dune::ParameterTree usedDefaultParams_;
};

} // end namespace Dumux

#endif
