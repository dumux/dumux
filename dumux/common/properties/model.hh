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
 * \ingroup Properties
 * \brief Defines a type tags and some fundamental properties for all models
 */
#ifndef DUMUX_MODEL_PROPERTIES_HH
#define DUMUX_MODEL_PROPERTIES_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/balanceequationopts.hh>
#include <dumux/io/defaultiofields.hh>

// Forward declaration
namespace Dune { class ParameterTree; }

namespace Dumux {
namespace Properties {

//! Type tag for numeric models.
namespace TTag {
struct ModelProperties {};
}

//! Set the default type of scalar values to double
template<class TypeTag>
struct Scalar<TypeTag, TTag::ModelProperties> { using type = double; };

//! Set the default vector with size number of equations to a field vector
template<class TypeTag>
struct NumEqVector<TypeTag, TTag::ModelProperties> { using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

//! Set the default primary variable vector to a vector of size of number of equations
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ModelProperties> { using type = GetPropType<TypeTag, Properties::NumEqVector>; };

//! do not specific any model-specific default parameters here
template<class TypeTag>
struct ModelDefaultParameters<TypeTag, TTag::ModelProperties>
{
    static void defaultParams(Dune::ParameterTree& tree, const std::string& group = "") { }
};

//! Set the default to an implementation throwing a NotImplemented error
template<class TypeTag>
struct IOFields<TypeTag, TTag::ModelProperties> { using type = DefaultIOFields; };

//! Set the default class for the balance equation options
template<class TypeTag>
struct BalanceEqOpts<TypeTag, TTag::ModelProperties> { using type = BalanceEquationOptions<TypeTag>; };

} // namespace Properties
} // namespace Dumux

#endif
