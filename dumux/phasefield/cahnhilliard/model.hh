// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CahnHilliardModel
 * \brief The Cahn-Hilliard equation
 *
 * The Cahn-Hilliard equation models phase separation and is given by
 \f[
 \frac{\partial c}{\partial t} = \nabla\cdot\left( M \nabla \mu \right)
 \f]
 \f[
 \mu = \frac{\partial f}{\partial c} - \gamma \nabla^2 c
 \f]
 * where \f$ c \f$ is the concentration, \f$ \mu \f$ is the chemical potential,
 * \f$ M \f$ is the mobility, \f$ \gamma \f$ is the surface tension (gradient energy
 * coefficient), and \f$ f(c) \f$ is a double-well free energy density. This model
 * splits the fourth-order equation into two coupled second-order equations with
 * \f$ c \f$ and \f$ \mu \f$ as primary variables, so that it can be discretized with
 * (piecewise linear or quadratic) CVFE schemes.
 */
#ifndef DUMUX_PHASEFIELD_CAHNHILLIARD_MODEL_HH
#define DUMUX_PHASEFIELD_CAHNHILLIARD_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>

#include "indices.hh"
#include "localresidual.hh"
#include "variables.hh"

namespace Dumux::Properties {

//! Type tag for the Cahn-Hilliard model
namespace TTag {
struct CahnHilliard { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! Use the local residual of the Cahn-Hilliard model
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::CahnHilliard>
{ using type = CahnHilliardLocalResidual<TypeTag>; };

//! The model traits of the Cahn-Hilliard model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::CahnHilliard>
{
    struct type
    {
        using Indices = CahnHilliardIndices;
        static constexpr int numEq() { return 2; }
    };
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::CahnHilliard>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = CahnHilliardVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
