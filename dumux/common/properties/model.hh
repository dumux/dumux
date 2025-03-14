// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
#include <dumux/discretization/defaultlocaloperator.hh>

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

//! Set the default primary variable vector to a vector of size of number of equations
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ModelProperties> { using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                                                         GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

//! Set the default to an implementation throwing a NotImplemented error
template<class TypeTag>
struct IOFields<TypeTag, TTag::ModelProperties> { using type = DefaultIOFields; };

//! Set the default class for the balance equation options
template<class TypeTag>
struct BalanceEqOpts<TypeTag, TTag::ModelProperties> { using type = BalanceEquationOptions<TypeTag>; };

template<class TypeTag>
class DeprecatedBaseLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    struct PropertyBaseLocalResidual {
        [[deprecated("BaseLocalResidual property is deprecated. Will be removed after release 3.10. Use DiscretizationDefaultLocalOperator.")]]
        PropertyBaseLocalResidual() = default;
        int dummy = 0;
    };
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
public:
    using ParentType::ParentType;
private:
    PropertyBaseLocalResidual deprecated_;
};

//! Deprecation helper for BaseLocalResidual
template<class TypeTag>
struct [[deprecated("BaseLocalResidual property is deprecated. Will be removed after release 3.10. Use DiscretizationDefaultLocalOperator.")]]
BaseLocalResidual<TypeTag, TTag::ModelProperties>
{
    using type [[deprecated("BaseLocalResidual property is deprecated. Will be removed after release 3.10. Use DiscretizationDefaultLocalOperator.")]]
        = DeprecatedBaseLocalResidual<TypeTag>;
};

} // namespace Properties
} // namespace Dumux

#endif
