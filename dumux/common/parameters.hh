/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief The infrastructure to retrieve run-time parameters from
 *        Dune::ParameterTrees with the defaul value taken from the
 *        property system.
 */
#ifndef DUMUX_PARAMETERS_HH
#define DUMUX_PARAMETERS_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/exceptions.hh>

#include <dune/common/parametertree.hh>

#include <sstream>
#include <list>
#include <unordered_map>

/*!
 * \brief Retrieve a runtime parameter which _does_ have a default value taken from
 *        the Dumux property system.
 *
 * If the macro is called with four argument the third argument is
 * group name which must be the prefix to the property name which
 * provides the default value for the parameter
 *
 * Examples:
 *
 * // -> retrieves scalar value UpwindWeight, default
 * // is taken from the property UpwindWeight
 * GET_PARAM(TypeTag, Scalar, UpwindWeight);
 *
 * // -> retrieves Boolean value Newton.WriteConvergence, default
 * // is taken from the property NewtonWriteConvergence
 * GET_PARAM(TypeTag, bool, Newton, WriteConvergence);
 */
#define GET_PARAM(TypeTag, ParamType, ParamNameOrGroupName, ...)         \
    Dumux::Parameters::get<TypeTag,                                     \
                           ParamType,                                   \
                           PTAG(ParamNameOrGroupName ## __VA_ARGS__)>   \
    (#ParamNameOrGroupName,                                              \
     Dumux::Parameters::getString_(#__VA_ARGS__))

/*!
 * \brief Retrieve a runtime parameter which _does not_ have a default value taken from
 *        the Dumux property system.
 *
 * If the macro is called with four argument the third argument is
 * group name
 *
 * Examples:
 *
 * // -> retrieves global integer value NumberOfCellsX
 * GET_RUNTIME_PARAM(TypeTag, int, NumberOfCellsX);
 *
 * // -> retrieves global integer value NumberOfCellsX which is
 * // located int the parameter group "Geometry"
 * GET_RUNTIME_PARAM(TypeTag, int, Geometry, NumberOfCellsX);
 */
#define GET_RUNTIME_PARAM(TypeTag, ParamType, ParamNameOrGroupName, ...) \
        Dumux::Parameters::getRuntime<TypeTag,                          \
                                      ParamType>                        \
     (#ParamNameOrGroupName,                                            \
      Dumux::Parameters::getString_(#__VA_ARGS__))

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
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;
    const Dune::ParameterTree &rt = Params::runTimeParams();

    // loop over all keys of the current tree
    const Dune::ParameterTree::KeyVector &keys = 
        tree.getValueKeys();
    for (int i = 0; i < keys.size(); ++i) {
        std::string canonicalName = prefix + keys[i];

        // check whether the key was accessed
        if (rt.hasKey(canonicalName))
            continue;
        unusedParams.push_back(canonicalName);
    }

    // loop over all subtrees
    const Dune::ParameterTree::KeyVector &subKeys = 
        tree.getSubKeys();
    for (int i = 0; i < subKeys.size(); ++i) {
        std::string newPrefix = prefix + subKeys[i] + ".";

        findUnusedKeys_<TypeTag>(unusedParams, 
                                 tree.sub(subKeys[i]),
                                 newPrefix);
    }
    
}

template <class TypeTag>
void print(std::ostream &os = std::cout)
{
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

    const Dune::ParameterTree &tree = Params::tree();
    const Dune::ParameterTree &rt = Params::runTimeParams();
    const Dune::ParameterTree &ct = Params::compileTimeParams();

    os << "###############################\n";
    os << "# Run-time parameters:\n";
    os << "###############################\n";
    rt.report(os);
    os << "###############################\n";
    os << "# Compile-time parameters:\n";
    os << "###############################\n";
    ct.report(os);

    std::list<std::string> unusedParams;
    findUnusedKeys_<TypeTag>(unusedParams, tree);

    if (unusedParams.size() > 0) {
        os << "###############################\n";
        os << "# UNUSED PARAMETERS:\n";
        os << "###############################\n";
        std::list<std::string>::const_iterator it = unusedParams.begin();
        for (; it != unusedParams.end(); ++it) {
            os << *it << " = \"" << tree.get(*it, "") << "\"\n";
        };
    }
};

const char *getString_(const char *foo = 0)
{ return foo; }

template <class TypeTag> 
class Param
{
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;
public:
    template <class ParamType, class PropTag>
    static const ParamType &get(const char *groupOrParamName, 
                                const char *paramNameOrNil = 0)
    {
#ifndef NDEBUG
        // make sure that the parameter is used consistently. since
        // this is potentially quite expensive, it is only done if
        // debugging code is not explicitly turned off.
        const char *paramName, *groupName;
        std::string propertyName;
        if (paramNameOrNil && strlen(paramNameOrNil) > 0) {
            groupName = groupOrParamName;
            paramName = paramNameOrNil;
            propertyName  = groupName;
            propertyName += paramName;
        }
        else {
            groupName = "";
            paramName = groupOrParamName;
            propertyName = paramName;
        }

        check_<ParamType>(propertyName, groupName, paramName);
#endif

        static const ParamType &value = retrieve_<ParamType, PropTag>(groupOrParamName, paramNameOrNil);
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

        check_<ParamType>(propertyName, groupName, paramName);
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
        };
    };

    template <class ParamType>
    static void check_(const std::string &propertyName, 
                       const char *groupName, 
                       const char *paramName)
    {
        const std::string &paramTypeName = 
            Dune::className<ParamType>();
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
    static const ParamType &retrieve_(const char *groupOrParamName, const char *paramNameOrNil = 0)
    {   
        const char *paramName, *groupName;
        if (paramNameOrNil && strlen(paramNameOrNil) > 0) {
            groupName = groupOrParamName;
            paramName = paramNameOrNil;
        }
        else {
            groupName = 0;
            paramName = groupOrParamName;
        }

        // prefix the parameter name by 'GroupName.'. E.g. 'Newton'
        // and 'WriteConvergence' becomes 'Newton.WriteConvergence'
        // with the default value specified by the
        // 'NewtonWriteConvergence' property. in an INI file this
        // would look like:
        //
        // [Newton]
        // WriteConvergence = true
        std::string canonicalName(paramName);
        if (groupName && strlen(groupName) > 0) {
            canonicalName.insert(0, ".");
            canonicalName.insert(0, groupName);
        }

        std::string modelParamGroup(GET_PROP(TypeTag, PTAG(ModelParameterGroup))::value());
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

        // retrieve actual parameter from the parameter tree
        ParamType defaultValue = GET_PROP_VALUE(TypeTag, PropTag);
        static ParamType value = Params::tree().template get<ParamType>(canonicalName, defaultValue);

        // remember whether the parameter was taken from the parameter
        // tree or the default from the property system was taken.
        Dune::ParameterTree &rt = Params::runTimeParams();
        Dune::ParameterTree &ct = Params::compileTimeParams();
        if (Params::tree().hasKey(canonicalName)) {
            rt[canonicalName] = Params::tree()[canonicalName];
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
            groupName = 0;
            paramName = groupOrParamName;
        }

        static std::string modelParamGroup(GET_PROP(TypeTag, PTAG(ModelParameterGroup))::value());

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
        if (groupName && groupName[0] != '\0') {
            canonicalName.append(groupName);
            canonicalName.push_back('.');
        }

        // append the name of the parameter
        canonicalName.append(paramName);

        // cache parameters using a hash_map (Dune::Parameter tree is slow!)
        typedef std::unordered_map<std::string, ParamType> ParamCache;
        static ParamCache paramCache;
        const typename ParamCache::iterator &it = paramCache.find(canonicalName);
        if (it != paramCache.end())
            return it->second;

        // retrieve actual parameter from the parameter tree
        if (!Params::tree().hasKey(canonicalName)) {
            DUNE_THROW(Dumux::ParameterException,
                       "Mandatory parameter '" << canonicalName
                       << "' was not specified");
        }

        // update the cache
        ParamType defaultValue;
        ParamType value = Params::tree().template get<ParamType>(canonicalName, defaultValue);
        paramCache[canonicalName] = value;

        // remember whether the parameter was taken from the parameter
        // tree or the default from the property system was taken.
        Dune::ParameterTree &rt = Params::runTimeParams();
        rt[canonicalName] = Params::tree()[canonicalName];

        return paramCache[canonicalName];
    }
};

template <class TypeTag, class ParamType, class PropTag>
const ParamType &get(const char *paramOrGroupName,
                     const char *paramNameOrNil = 0)
{
    return Param<TypeTag>::template get<ParamType, PropTag>(paramOrGroupName,
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
