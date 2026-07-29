// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Fluid system for the classic (purely diffusive) Henry problem, parameterized
 *        to match "Test Case 1" of Fahs et al. (2016, WRR,
 *        doi:10.1002/2016WR019288) -- their semianalytical validation case, itself the
 *        original Henry (1964) configuration with an exaggerated molecular diffusion
 *        coefficient.
 */
#ifndef DUMUX_TEST_HENRY_FAHS_FLUIDSYSTEM_HH
#define DUMUX_TEST_HENRY_FAHS_FLUIDSYSTEM_HH

#include <string>
#include <limits>

#include <dumux/io/name.hh>
#include <dumux/material/fluidsystems/base.hh>

namespace Dumux::FluidSystems {

/*!
 * \ingroup OnePNCTests
 * \brief Single-phase, two-component fluid system matching Fahs et al. (2016) Test
 *        Case 1 (Table 2): rho0 = 1000 kg/m^3, rho1 = 1025 kg/m^3, Dm = 18.86e-6 m^2/s,
 *        no dispersion (alphaL = alphaT = 0, i.e. purely diffusive -- matching that
 *        this model has no dispersion tensor at all).
 *
 * The salt mass fraction is scaled so that X = 0.035 at the seaward boundary; only the
 * ratio X/0.035 (the paper's dimensionless concentration c) is physically meaningful
 * for comparison against the paper's results.
 */
template<class Scalar>
class HenryFahsFluid : public Base<Scalar, HenryFahsFluid<Scalar>>
{
    using ThisType = HenryFahsFluid<Scalar>;

public:
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 2;

    static constexpr int phase0Idx = 0;
    static constexpr int liquidPhaseIdx = phase0Idx;

    static constexpr int solventIdx = 0;
    static constexpr int soluteIdx = 1;
    static constexpr int comp0Idx = solventIdx;
    static constexpr int comp1Idx = soluteIdx;

    //! Freshwater density [kg/m^3]
    static constexpr Scalar rhoFresh()
    { return 1000.0; }

    //! Seawater density [kg/m^3] (Fahs et al. 2016, Table 2: rho1 = 1.025e3 kg/m^3)
    static constexpr Scalar rhoSeawater()
    { return 1025.0; }

    //! Seawater salt mass fraction [-]: an arbitrary (but fixed) scale -- only
    //! X/seawaterMassFraction() (i.e. the dimensionless concentration c) is compared
    //! against the paper's results.
    static constexpr Scalar seawaterMassFraction()
    { return 0.035; }

    static std::string phaseName(int phaseIdx = liquidPhaseIdx)
    { return IOName::liquidPhase(); }

    static constexpr bool isMiscible()
    { return false; }

    static constexpr bool isGas(int phaseIdx = liquidPhaseIdx)
    { return false; }

    static constexpr bool isIdealMixture(int phaseIdx = liquidPhaseIdx)
    { return true; }

    static constexpr bool isCompressible(int phaseIdx = liquidPhaseIdx)
    { return false; } // density does not depend on pressure

    static constexpr bool isIdealGas(int phaseIdx = liquidPhaseIdx)
    { return false; }

    static std::string componentName(int compIdx)
    { return compIdx == solventIdx ? "solvent" : "solute"; }

    static Scalar molarMass(int compIdx)
    { return compIdx == solventIdx ? 18.015e-3 : 58.44e-3; } // H2O, NaCl [kg/mol]

    static void init() {}

    using Base<Scalar, ThisType>::density;
    /*!
     * \brief rho(X) = rhoFresh + (rhoSeawater - rhoFresh)/Xsw * X
     */
    template<class FluidState>
    static Scalar density(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        static constexpr Scalar densitySlope = (rhoSeawater() - rhoFresh())/seawaterMassFraction();
        return rhoFresh() + densitySlope*fluidState.massFraction(liquidPhaseIdx, soluteIdx);
    }

    using Base<Scalar, ThisType>::molarDensity;
    //! \copydoc Base<Scalar,ThisType>::molarDensity(const FluidState&,int)
    template<class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    { return density(fluidState, phaseIdx)/molarMass(solventIdx); }

    using Base<Scalar, ThisType>::viscosity;
    //! Constant viscosity (Fahs et al. 2016, Table 2: mu = 1e-3 kg/(m s))
    template<class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    { return 1e-3; }

    using Base<Scalar, ThisType>::fugacityCoefficient;
    template<class FluidState>
    static Scalar fugacityCoefficient(const FluidState& fluidState, int phaseIdx, int compIdx)
    { return phaseIdx == compIdx ? 1.0 : std::numeric_limits<Scalar>::infinity(); }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    //! Molecular diffusion coefficient Dm = 18.86e-6 m^2/s (Fahs et al. 2016, Table 2,
    //! Test Case 1 -- the exaggerated value matching the original Henry [1964] problem)
    template<class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx, int compIIdx, int compJIdx)
    { return 18.86e-6; }
};

} // end namespace Dumux::FluidSystems

#endif
