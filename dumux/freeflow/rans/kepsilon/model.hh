// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief The two-equation, high-Reynolds-number wall-function k-epsilon RANS turbulence
 *        closure, fused as two extra transported equations (turbulent kinetic energy k,
 *        dissipation rate epsilon) onto the single-phase Navier-Stokes mass balance.
 *
 * See turbulenceequations.md \S8 for the physics, and whatisimplemented.md/
 * proposedimplementation.md for why these equations live on the *mass* sub-model, following
 * the same strategy already validated for the one-equation and k-omega models.
 */
#ifndef DUMUX_RANS_KEPSILON_MASS_MODEL_HH
#define DUMUX_RANS_KEPSILON_MASS_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "localresidual.hh"
#include "iofields.hh"
#include "advectiveflux.hh"

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Traits for the k-epsilon turbulence closure fused onto the single-phase Navier-Stokes
 *        mass model.
 */
struct KEpsilonMassModelTraits : public NavierStokesMassOnePModelTraits
{
    //! One mass balance equation, plus two turbulence transport equations (k, epsilon).
    static constexpr int numEq() { return 3; }

    using Indices = KEpsilonMassIndices;
};

namespace Properties {

namespace TTag {
//! The type tag for the single-phase Navier-Stokes mass model fused with the two-equation
//! k-epsilon turbulence closure.
struct NavierStokesMassOneKEpsilon { using InheritsFrom = std::tuple<NavierStokesMassOneP>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOneKEpsilon>
{ using type = KEpsilonMassModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassOneKEpsilon>
{ using type = KEpsilonMassLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOneKEpsilon>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMassOnePVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = KEpsilonMassVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOneKEpsilon> { using type = KEpsilonMassIOFields; };

} // end namespace Properties
} // end namespace Dumux

#endif
