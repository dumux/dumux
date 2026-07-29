// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Artificial fluid system for the lock-exchange test.
 */
#ifndef DUMUX_TEST_LOCKEXCHANGE_FLUIDSYSTEM_HH
#define DUMUX_TEST_LOCKEXCHANGE_FLUIDSYSTEM_HH

#include <string>
#include <limits>

#include <dumux/io/name.hh>
#include <dumux/material/fluidsystems/base.hh>

namespace Dumux::FluidSystems {

/*!
 * \ingroup OnePNCTests
 * \brief Artificial single-phase, two-component fluid system for the lock-exchange
 *        test: density is an exactly linear function of the tracer mass fraction,
 *        rho = rho0*(1 + C). Deliberately not based on any real component data (unlike
 *        e.g. FluidSystems::Brine) so that the density variance driving the lock
 *        exchange can be dialed to a small, known value (a couple of percent) that
 *        stays within the Boussinesq approximation's validity range -- the point of
 *        this test is to check the approximation where it is expected to hold,
 *        as opposed to the saltwaterintrusion test's realistic (and much larger, up
 *        to 3.5%) salinity range, where it does not hold well.
 *
 * \tparam Scalar scalar type
 * \tparam variableDensity if true, density varies with the tracer mass fraction (the
 *         "full" model); if false, density is the constant reference density (the
 *         Boussinesq model, where the tracer's effect on density re-enters only
 *         through Problem::buoyantDensity() in BoussinesqCVFEDarcyLaw)
 */
template<class Scalar, bool variableDensity>
class LockExchangeFluid : public Base<Scalar, LockExchangeFluid<Scalar, variableDensity>>
{
    using ThisType = LockExchangeFluid<Scalar, variableDensity>;

public:
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 2;

    static constexpr int phase0Idx = 0;
    static constexpr int liquidPhaseIdx = phase0Idx;

    static constexpr int solventIdx = 0;
    static constexpr int soluteIdx = 1;
    static constexpr int comp0Idx = solventIdx;
    static constexpr int comp1Idx = soluteIdx;

    //! The constant reference density [kg/m^3]: the density at zero tracer concentration
    static constexpr Scalar referenceDensity()
    { return 1000.0; }

    static std::string phaseName(int phaseIdx = liquidPhaseIdx)
    { return IOName::liquidPhase(); }

    static constexpr bool isMiscible()
    { return false; }

    static constexpr bool isGas(int phaseIdx = liquidPhaseIdx)
    { return false; }

    static constexpr bool isIdealMixture(int phaseIdx = liquidPhaseIdx)
    { return true; }

    static constexpr bool isCompressible(int phaseIdx = liquidPhaseIdx)
    { return false; } // density does not depend on pressure in either variant

    static constexpr bool isIdealGas(int phaseIdx = liquidPhaseIdx)
    { return false; }

    static std::string componentName(int compIdx)
    { return compIdx == solventIdx ? "solvent" : "solute"; }

    static Scalar molarMass(int compIdx)
    { return 18e-3; } // both components given the same (arbitrary) molar mass; only mass fractions matter here

    static void init() {}

    using Base<Scalar, ThisType>::density;
    /*!
     * \brief rho = rho0*(1+C) if variableDensity, else the constant reference rho0.
     */
    template<class FluidState>
    static Scalar density(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        if constexpr (variableDensity)
            return referenceDensity()*(1.0 + fluidState.massFraction(liquidPhaseIdx, soluteIdx));
        else
            return referenceDensity();
    }

    using Base<Scalar, ThisType>::molarDensity;
    //! \copydoc Base<Scalar,ThisType>::molarDensity(const FluidState&,int)
    template<class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    { return density(fluidState, phaseIdx)/molarMass(solventIdx); }

    using Base<Scalar, ThisType>::viscosity;
    template<class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    { return 1e-3; } // water-like, constant

    using Base<Scalar, ThisType>::fugacityCoefficient;
    template<class FluidState>
    static Scalar fugacityCoefficient(const FluidState& fluidState, int phaseIdx, int compIdx)
    { return phaseIdx == compIdx ? 1.0 : std::numeric_limits<Scalar>::infinity(); }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    template<class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx, int compIIdx, int compJIdx)
    { return 1e-9; }
};

} // end namespace Dumux::FluidSystems

#endif
