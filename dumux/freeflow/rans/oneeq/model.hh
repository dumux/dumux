// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief The one-equation (Spalart-Allmaras) RANS turbulence closure, fused as a second
 *        transported equation (the working viscosity ν̃) onto the single-phase Navier-Stokes
 *        mass balance.
 *
 * See turbulenceequations.md \S7 for the physics, and whatisimplemented.md/
 * proposedimplementation.md for why this equation lives on the *mass* sub-model (mirroring
 * how releases/3.10 had ν̃ share the same monolithic primary-variable block as pressure and
 * velocity) rather than on a dedicated third coupled sub-domain.
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_MODEL_HH
#define DUMUX_RANS_ONEEQ_MASS_MODEL_HH

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
 * \brief Traits for the one-equation (Spalart-Allmaras) turbulence closure fused onto the
 *        single-phase Navier-Stokes mass model.
 */
struct OneEqMassModelTraits : public NavierStokesMassOnePModelTraits
{
    //! One mass balance equation, plus one turbulence transport equation (ν̃).
    static constexpr int numEq() { return 2; }

    using Indices = OneEqMassIndices;
};

namespace Properties {

namespace TTag {
//! The type tag for the single-phase Navier-Stokes mass model fused with the one-equation
//! (Spalart-Allmaras) turbulence closure.
struct NavierStokesMassOnePOneEq { using InheritsFrom = std::tuple<NavierStokesMassOneP>; };
} // end namespace TTag

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NavierStokesMassOnePOneEq>
{ using type = OneEqMassModelTraits; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NavierStokesMassOnePOneEq>
{ using type = OneEqMassLocalResidual<TypeTag>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::NavierStokesMassOnePOneEq>
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
    using type = OneEqMassVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::NavierStokesMassOnePOneEq> { using type = OneEqMassIOFields; };

} // end namespace Properties
} // end namespace Dumux

#endif
