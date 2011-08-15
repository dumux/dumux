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

#include "propertysystem.hh"

#include <dune/common/parametertree.hh>

#include <sstream>
#include <list>

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
#define GET_PARAM(TypeTag, ParamType, ParamNameOrGroupName, ...) \
    Dumux::Parameters::Param<TypeTag, PTAG(ParamNameOrGroupName ## __VA_ARGS__), ParamType> \
    :: get(#ParamNameOrGroupName, Dumux::Parameters::getString_(#__VA_ARGS__))

// retrieve a parameter which does _not_ have a default value taken
// from the dumux property system
//#define GET_RUNTIME_PARAM(TypeTag, ParamType, ParamName) 
//    Dumux::Parameters::RunTimeParam<TypeTag, ParamType>::get(#ParamName)

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(ParameterTree);
NEW_PROP_TAG(ModelParameterGroup);
} // namespace Properties

namespace Parameters {
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
    int n = 0;
    const Dune::ParameterTree::KeyVector &keys = 
        Params::tree().getValueKeys();
    for (int i = 0; i < keys.size(); ++i) {
        // check wheter the key was accessed
        if (rt.hasKey(keys[i]))
            continue;
        ++n;
        unusedParams.push_back(keys[i]);
    }

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

template <class TypeTag, 
          class PropTag,
          class ParamType> 
class Param
{
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;
public:
    static const ParamType &get(const char *groupOrParamName, const char *paramNameOrNil = 0)
    {
#ifndef NDEBUG
        const char *paramName, *groupName;
        if (paramNameOrNil) {
            groupName = groupOrParamName;
            paramName = paramNameOrNil;
        }
        else {
            groupName = "";
            paramName = groupOrParamName;
        }
        // make sure that a property is only the default for a single
        // parameter, i.e. that a property is not used for multiple
        // parameters
        static std::string fixedGroupName = groupName;
        if (fixedGroupName != groupName)
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "Conflicting prefixes for parameter '" << paramName
                       << "': '" << fixedGroupName
                       << "' and '" << groupName << "'.");
        }
#endif

        static const ParamType &value = retrieve_(groupOrParamName, paramNameOrNil);
        return value;
    }

private:
    static const ParamType &retrieve_(const char *groupOrParamName, const char *paramNameOrNil = 0)
    {   
        const char *paramName, *groupName;
        if (paramNameOrNil) {
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
        std::string finalName(paramName);
        if (groupName) {
            finalName.insert(0, ".");
            finalName.insert(0, groupName);
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
            finalName.insert(0, ".");
            finalName.insert(0, modelParamGroup);
        }
        
        // retrieve actual parameter from the parameter tree
        ParamType defaultValue = GET_PROP_VALUE(TypeTag, PropTag);
        static ParamType value = Params::tree().template get<ParamType>(finalName, defaultValue);

        // remember whether the parameter was taken from the parameter
        // tree or the default from the property system was taken.
        Dune::ParameterTree &rt = Params::runTimeParams();
        Dune::ParameterTree &ct = Params::compileTimeParams();
        if (Params::tree().hasKey(finalName)) {
            rt[finalName] = Params::tree()[finalName];
        }
        else {
            std::string s;
            std::ostringstream oss(s);
            oss << defaultValue;
            ct[finalName] = oss.str();
        }
        return value;
    }
};

template <class TypeTag, 
          class ParamType> 
class RunTimeParam
{
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

public:
    static const ParamType &get(const char *name)
    {
        static const ParamType &value = retrieve_(name);
        return value;
    }

private:
    static const ParamType &retrieve_(const char *name)
    {   
        static ParamType value = Params::tree().template get<ParamType>(name);

        Dune::ParameterTree &rt = Params::tree().sub("RunTimeParams");
        if (Params::tree().hasKey(name)) {
            rt[name] = Params::tree()[name];
        }

        return value;
    }
};

} // namespace Parameters

} // namespace Dumux


#endif
