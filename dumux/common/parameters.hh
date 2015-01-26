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

#include <dune/common/parametertree.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/exceptions.hh>

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
#define GET_PARAM(TypeTag, ParamType, ParamName)                        \
    ::Dumux::Parameters::get<TypeTag,                                   \
                           ParamType,                                   \
                           PTAG_(ParamName)>(#ParamName, #ParamName)

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
    ::Dumux::Parameters::get<TypeTag,                                   \
                           ParamType,                                   \
                           PTAG_(GroupName##ParamName)>(#GroupName#ParamName, #GroupName, #ParamName)

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
    ::Dumux::Parameters::getRuntime<TypeTag, ParamType>(#ParamName)

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
    ::Dumux::Parameters::getRuntime<TypeTag, ParamType>(#GroupName, #ParamName)

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(ParameterTree);
NEW_PROP_TAG(ModelParameterGroup);
} // namespace Properties

namespace Parameters {

template <class TypeTag>
void findUnusedKeys_(std::list<std::string> &unusedParams,
                     const Dune::ParameterTree &tree,
                     const std::string prefix="")
{
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;
    const Dune::ParameterTree &rt = Params::runTimeParams();
    const Dune::ParameterTree &drt = Params::deprecatedRunTimeParams();

    // loop over all keys of the current tree
    const Dune::ParameterTree::KeyVector &keys =
        tree.getValueKeys();
    for (unsigned int i = 0; i < keys.size(); ++i) {
        std::string canonicalName = prefix + keys[i];

        // store keys which were not accessed
        if (!rt.hasKey(canonicalName) && !drt.hasKey(canonicalName))
        {
            unusedParams.push_back(canonicalName);
        }
    }

    // loop over all subtrees
    const Dune::ParameterTree::KeyVector &subKeys =
        tree.getSubKeys();
    for (unsigned int i = 0; i < subKeys.size(); ++i) {
        std::string newPrefix = prefix + subKeys[i] + ".";

        findUnusedKeys_<TypeTag>(unusedParams,
                                 tree.sub(subKeys[i]),
                                 newPrefix);
    }

}

template <class TypeTag>
bool hasDeprecatedKeys_(const Dune::ParameterTree &tree)
{
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;
    const Dune::ParameterTree &drt = Params::deprecatedRunTimeParams();

    // loop over all keys of the current tree
    const Dune::ParameterTree::KeyVector &keys =
        tree.getValueKeys();
    for (unsigned int i = 0; i < keys.size(); ++i) {
        std::string canonicalName = keys[i];

        // check whether the key was accessed
        if (drt.hasKey(canonicalName))
            return true;
    }
    return false;
}

/*!
 * \ingroup Parameter
 * \brief Print the run- and compile-time parameters.
 */
template <class TypeTag>
void print(std::ostream &os = std::cout)
{
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;

    const Dune::ParameterTree &tree = Params::tree();
    const Dune::ParameterTree &rt = Params::runTimeParams();
    const Dune::ParameterTree &ct = Params::compileTimeParams();
    const Dune::ParameterTree &drt = Params::deprecatedRunTimeParams();
    const Dune::ParameterTree &unrt = Params::unusedNewRunTimeParams();

    os << "# Run-time specified parameters:" << std::endl;
    rt.report(os);

    if (hasDeprecatedKeys_<TypeTag>(tree)) 
    {
        os << "# DEPRECATED Run-time specified parameters:" << std::endl;
        drt.report(os);
        os << "# Replace by:" << std::endl;
        unrt.report(os);
    }

    os << "# Compile-time specified parameters:" << std::endl;
    ct.report(os);

    std::list<std::string> unusedParams;
    findUnusedKeys_<TypeTag>(unusedParams, tree);

    if (unusedParams.size() > 0)
    {
        os << "# UNUSED PARAMETERS:" << std::endl;
        for (auto it = unusedParams.begin(); it != unusedParams.end(); ++it)
        {
            os << *it << " = \"" << tree.get(*it, "") << "\"" << std::endl;
        }
    }
}

const char *getString_(const char *foo = 0)
{ return foo; }

template <class TypeTag>
class Param
{
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;
public:
    template <class ParamType, class PropTag>
    static const ParamType &get(const char *propertyName,
                                const char *groupOrParamName,
                                const char *paramNameOrNil = 0)
    {
        static const ParamType &value = retrieve_<ParamType, PropTag>(propertyName, groupOrParamName, paramNameOrNil);
        return value;
    }

    template <class ParamType>
    static const ParamType &getRuntime(const char *groupOrParamName,
                                       const char *paramNameOrNil = 0)
    {
#ifndef NDEBUG
        // make sure that the parameter is used consistently. since
        // this is potentially quite expensive, it is only done if
        // debugging code is not explicitly turned off.
        const char *paramName, *groupName;
        static const std::string propertyName("");
        if (paramNameOrNil && strlen(paramNameOrNil) > 0) {
            groupName = groupOrParamName;
            paramName = paramNameOrNil;
        }
        else {
            groupName = "";
            paramName = groupOrParamName;
        }

        check_(Dune::className<ParamType>(), propertyName, groupName, paramName);
#endif

        return retrieveRuntime_<ParamType>(groupOrParamName, paramNameOrNil);
    }

private:
    struct Blubb {
        std::string propertyName;
        std::string paramTypeName;
        std::string groupName;

        Blubb &operator=(const Blubb &b)
        {
            propertyName = b.propertyName;
            paramTypeName = b.paramTypeName;
            groupName = b.groupName;
            return *this;
        }
    };

    static void check_(const std::string &paramTypeName,
                       const std::string &propertyName,
                       const char *groupName,
                       const char *paramName)
    {
        typedef std::unordered_map<std::string, Blubb> StaticData;
        static StaticData staticData;

        typename StaticData::iterator it = staticData.find(paramName);
        Blubb *b;
        if (it == staticData.end())
        {
            Blubb a;
            a.propertyName = propertyName;
            a.paramTypeName = paramTypeName;
            a.groupName = groupName;
            staticData[paramName] = a;
            b = &staticData[paramName];
        }
        else
            b = &(it->second);

        if (b->groupName != groupName) {
            DUNE_THROW(Dune::InvalidStateException,
                       "GET_*_PARAM for parameter '" << paramName
                       << "' called for at least two different groups ('"
                       << b->groupName << "' and '" << groupName << "')");
        }

        if (b->propertyName != propertyName) {
            DUNE_THROW(Dune::InvalidStateException,
                       "GET_*_PARAM for parameter '" << paramName
                       << "' called for at least two different properties ('"
                       << b->propertyName << "' and '" << propertyName << "')");
        }

        if (b->paramTypeName != paramTypeName) {
            DUNE_THROW(Dune::InvalidStateException,
                       "GET_*_PARAM for parameter '" << paramName << "' in group '"
                       << groupName << "' called with at least two different types ("
                       << b->paramTypeName << " and " << paramTypeName << ")");
        }
    }

    template <class ParamType, class PropTag>
    static const ParamType &retrieve_(const char *propertyName,
                                      const char *groupOrParamName,
                                      const char *paramNameOrNil = 0)
    {
        const char *paramName, *groupName;
        if (paramNameOrNil && strlen(paramNameOrNil) > 0) {
            groupName = groupOrParamName;
            paramName = paramNameOrNil;
        }
        else {
            groupName = "";
            paramName = groupOrParamName;
        }

#ifndef NDEBUG
        // make sure that the parameter is used consistently. since
        // this is potentially quite expensive, it is only done if
        // debugging code is not explicitly turned off.
        check_(Dune::className<ParamType>(), propertyName, groupName, paramName);
#endif

        // prefix the parameter name by 'GroupName.'. E.g. 'Newton'
        // and 'WriteConvergence' becomes 'Newton.WriteConvergence'
        // with the default value specified by the
        // 'NewtonWriteConvergence' property. in an INI file this
        // would look like:
        //
        // [Newton]
        // WriteConvergence = true
        std::string canonicalName(paramName);
        if (strlen(groupName) > 0) {
            canonicalName.insert(0, ".");
            canonicalName.insert(0, groupName);
        }

        std::string modelParamGroup(GET_PROP_VALUE(TypeTag, ModelParameterGroup));
        // prefix the parameter with the parameter group of the
        // model. this allows things like sub-model specific parameters like
        //
        // [Stokes.Newton]
        // WriteConvergence = false
        // [Darcy.Newton]
        // WriteConvergence = true
        if (modelParamGroup.size()) {
            canonicalName.insert(0, ".");
            canonicalName.insert(0, modelParamGroup);
        }

        static ParamType value;
        // retrieve actual parameter from the parameter tree
        ParamType defaultValue = GET_PROP_VALUE_(TypeTag, PropTag);
        if (!Params::tree().hasKey(canonicalName) && Params::tree().hasKey(paramName))//functionality to catch deprecated params
        {
            value = Params::tree().template get<ParamType>(paramName, defaultValue);
//            std::cout<<"\nWarning: Using the parameter: "<<paramName<<" without group name: "<<groupName<<" is deprecated!"<<"\n\n";
        }
        else
            value = Params::tree().template get<ParamType>(canonicalName, defaultValue);

        // remember whether the parameter was taken from the parameter
        // tree or the default from the property system was taken.
        Dune::ParameterTree &rt = Params::runTimeParams();
        Dune::ParameterTree &ct = Params::compileTimeParams();
        Dune::ParameterTree &drt = Params::deprecatedRunTimeParams();
        Dune::ParameterTree &unrt = Params::unusedNewRunTimeParams();
        if (Params::tree().hasKey(canonicalName)) {
            rt[canonicalName] = Params::tree()[canonicalName];
        }
        else if (Params::tree().hasKey(paramName))//functionality to catch deprecated params
        {
            drt[paramName] = Params::tree()[paramName];
            unrt[canonicalName] = Params::tree()[paramName];
        }
        else {
            std::string s;
            std::ostringstream oss(s);
            oss << defaultValue;
            ct[canonicalName] = oss.str();
        }
        return value;
    }

    template <class ParamType>
    static const ParamType &retrieveRuntime_(const char *groupOrParamName, const char *paramNameOrNil = 0)
    {
        const char *paramName, *groupName;
        if (paramNameOrNil && paramNameOrNil[0] != '\0') {
            groupName = groupOrParamName;
            paramName = paramNameOrNil;
        }
        else {
            groupName = "";
            paramName = groupOrParamName;
        }

        static std::string modelParamGroup(GET_PROP(TypeTag, ModelParameterGroup)::value);

        std::string canonicalName(modelParamGroup);

        // prefix the parameter with the parameter group of the
        // model. this allows things like sub-model specific parameters like
        //
        // [Stokes.Newton]
        // WriteConvergence = false
        // [Darcy.Newton]
        // WriteConvergence = true
        if (modelParamGroup.size()) {
            canonicalName.push_back('.');
        }

        // prefix the parameter name by 'GroupName.'. E.g. 'Newton'
        // and 'WriteConvergence' becomes 'Newton.WriteConvergence'
        // with the default value specified by the
        // 'NewtonWriteConvergence' property. in an INI file this
        // would look like:
        //
        // [Newton]
        // WriteConvergence = true
        if (strlen(groupName) > 0) {
            canonicalName.append(groupName);
            canonicalName.push_back('.');
        }

        // append the name of the parameter
        canonicalName.append(paramName);

        // cache parameters using a hash_map (Dune::Parameter tree is slow!)
        typedef std::unordered_map<std::string, ParamType> ParamCache;
        static ParamCache paramCache;
        typename ParamCache::iterator it = paramCache.find(canonicalName);
        if (it != paramCache.end())
            return it->second;

        it = paramCache.find(paramName);
        if (it != paramCache.end())
                    return it->second;

        // retrieve actual parameter from the parameter tree
        if (!Params::tree().hasKey(canonicalName) && !Params::tree().hasKey(paramName)) {
            print<TypeTag>();
            DUNE_THROW(::Dumux::ParameterException,
                       "Mandatory parameter '" << canonicalName
                       << "' was not specified");
        }

        // update the cache
        ParamType value;
        if (!Params::tree().hasKey(canonicalName) && Params::tree().hasKey(paramName))//functionality to catch deprecated params
        {
            value = Params::tree().template get<ParamType>(paramName);
            paramCache[paramName] = value;

            // remember whether the parameter was taken from the parameter
            // tree or the default from the property system was taken.
            Dune::ParameterTree &drt = Params::deprecatedRunTimeParams();
            Dune::ParameterTree &unrt = Params::unusedNewRunTimeParams();

            drt[paramName] = Params::tree()[paramName];
            unrt[canonicalName] = Params::tree()[paramName];
            return paramCache[paramName];
        }
        else
        {
            value = Params::tree().template get<ParamType>(canonicalName);
            paramCache[canonicalName] = value;

            // remember whether the parameter was taken from the parameter
            // tree or the default from the property system was taken.
            Dune::ParameterTree &rt = Params::runTimeParams();

            rt[canonicalName] = Params::tree()[canonicalName];
            return paramCache[canonicalName];
        }
    }
};

template <class TypeTag, class ParamType, class PropTag>
const ParamType &get(const char *propertyName,
                     const char *paramOrGroupName,
                     const char *paramNameOrNil = 0)
{
    return Param<TypeTag>::template get<ParamType, PropTag>(propertyName,
                                                            paramOrGroupName,
                                                            paramNameOrNil);
}

template <class TypeTag, class ParamType>
const ParamType &getRuntime(const char *paramOrGroupName,
                            const char *paramNameOrNil = 0)
{
    return Param<TypeTag>::template getRuntime<ParamType>(paramOrGroupName,
                                                          paramNameOrNil);
}

} // namespace Parameters

} // namespace Dumux


#endif
