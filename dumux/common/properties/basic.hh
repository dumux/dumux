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
 * \ingroup Properties
 * \file
 *
 * \brief Defines a type tags and some fundamental properties for
 *        all models
 */
#ifndef DUMUX_BASIC_PROPERTIES_HH
#define DUMUX_BASIC_PROPERTIES_HH

#include <dumux/common/balanceequationopts.hh>
#include <dumux/common/properties.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/defaultvtkoutputfields.hh>


namespace Dumux
{
namespace Properties
{
//! Type tag for numeric models.
NEW_TYPE_TAG(BasicProperties);

//! Set the default type of scalar values to double
SET_TYPE_PROP(BasicProperties, Scalar, double);

//! Set the default number of equations to one
SET_INT_PROP(BasicProperties, NumEq, 1);

//! use the global group as default for the model's parameter group
SET_STRING_PROP(BasicProperties, ModelParameterGroup, "");

//! do not specific any model-specific default parameters here
SET_PROP(BasicProperties, ModelDefaultParameters)
{
    static void defaultParams(Dune::ParameterTree& tree, const std::string& group = "") { }
};

//! Set the default to a function throwing a NotImplemented error
SET_TYPE_PROP(BasicProperties, VtkOutputFields, DefaultVtkOutputFields);

//! Set the default class for the balance equation options
SET_TYPE_PROP(BasicProperties, BalanceEqOpts, BalanceEquationOptions<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
