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

#define GET_PARAM(TypeTag, ParamType, ParamName)                         \
    Dumux::Parameters::Param<TypeTag, PTAG(ParamName), ParamType>::get(#ParamName)

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(ParameterTree);
} // namespace Properties

namespace Parameters {

template <class TypeTag, 
          class PropTag,
          class ParamType> 
class Param
{
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;
public:
    static const ParamType &get(const char *name)
    {
        static ParamType defaultValue = GET_PROP_VALUE(TypeTag, PropTag);
        static ParamType value = Params::tree().get(name, defaultValue);
        return value;
    }
};



} // namespace Parameters

} // namespace Dumux


#endif
