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

#include <dune/common/parametertree.hh>

#include <sstream>

#include "propertysystem.hh"

// retrieve a parameter which _does_ have a default value taken from
// the dumux property system. the optional last argument is the group
// name which must be the prefix to the property name which provides
// the default value for the parameter.
#define GET_PARAM(TypeTag, ParamType, ParamName, ...)      \
    Dumux::Parameters::Param<TypeTag, PTAG(__VA_ARGS__ ## ParamName), ParamType> \
    :: get(#ParamName, Dumux::Parameters::getGroupName_(#__VA_ARGS__))

// retrieve a parameter which does _not_ have a default value taken
// from the dumux property system
//#define GET_RUNTIME_PARAM(TypeTag, ParamType, ParamName)              \
//    Dumux::Parameters::RunTimeParam<TypeTag, ParamType>::get(#ParamName)

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(ParameterTree);
NEW_PROP_TAG(ModelParameterGroup);
} // namespace Properties

namespace Parameters {

const char *getGroupName_(const char *group = "")
{ return group; }

template <class TypeTag, 
          class PropTag,
          class ParamType> 
class Param
{
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;
public:
    static const ParamType &get(const char *paramName, const char *groupName = "")
    {
#ifndef NDEBUG
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

        static const ParamType &value = retrieve_(paramName, groupName);
        return value;
    }

private:
    static const ParamType &retrieve_(const char *paramName, const char *groupName)
    {   
        // prefix the parameter name by 'GroupName.'. E.g. 'Newton'
        // and 'WriteConvergence' becomes 'Newton.WriteConvergence'
        // with the default value specified by the
        // 'NewtonWriteConvergence' property. in an INI file this
        // would look like:
        //
        // [Newton]
        // WriteConvergence = true
        std::string finalName(paramName);
        int groupNameLen = strlen(groupName);
        if (groupNameLen) {
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
        Dune::ParameterTree &rt = Params::tree().sub("RunTimeParams");
        Dune::ParameterTree &ct = Params::tree().sub("DefaultParams");
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
