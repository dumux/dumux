// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \brief Liquid with constant properties and the matching fluid state, both without pressure work.
 *
 * The analytical solution of the convection test assumes \f$ h = c T \f$, i.e. it neglects the
 * pressure work. Here, the enthalpy equals the internal energy, in the fluid and in the fluid state.
 */

#ifndef DUMUX_1PNI_CONVECTION_CONSTANT_FLUID_NO_PRESSURE_WORK_HH
#define DUMUX_1PNI_CONVECTION_CONSTANT_FLUID_NO_PRESSURE_WORK_HH

#include <dumux/material/components/constant.hh>
#include <dumux/material/components/componenttraits.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux {

/*!
 * \brief Liquid with constant properties whose enthalpy is the internal energy.
 *
 * The enthalpy of Components::Constant is \f$ h = u + p/\rho \f$.
 */
template<int id, class Scalar>
class ConstantLiquidNoPressureWork : public Components::Constant<id, Scalar>
{
    using ParentType = Components::Constant<id, Scalar>;

public:
    //! Specific enthalpy of the liquid, without the pressure work
    static Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { return ParentType::liquidInternalEnergy(temperature, pressure); }
};

//! The traits of the derived component, as the default traits expect a CRTP base class
template<int id, class ScalarT>
struct ComponentTraits<ConstantLiquidNoPressureWork<id, ScalarT>>
{
    using Scalar = ScalarT;
    static constexpr bool hasSolidState = false;
    static constexpr bool hasLiquidState = true;
    static constexpr bool hasGasState = false;
    static constexpr bool isIon = false;
};

/*!
 * \brief Fluid state whose internal energy is the enthalpy.
 *
 * The internal energy of ImmiscibleFluidState is \f$ u = h - p/\rho \f$.
 */
template<class Scalar, class FluidSystem>
class ImmiscibleFluidStateNoPressureWork : public ImmiscibleFluidState<Scalar, FluidSystem>
{
    using ParentType = ImmiscibleFluidState<Scalar, FluidSystem>;

public:
    using ParentType::ParentType;

    //! The specific internal energy, which is the enthalpy
    Scalar internalEnergy(int phaseIdx) const
    { return this->enthalpy(phaseIdx); }
};

} // end namespace Dumux

#endif
