// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Fluid system for the Henry problem, parameterized to match the test cases of
 *        Fahs et al. (2016, WRR, doi:10.1002/2016WR019288). Density, viscosity and
 *        component naming are fixed (identical across all three of the paper's test
 *        cases, Table 2); only the molecular diffusion coefficient is a runtime
 *        parameter (Problem.MolecularDiffusion), since it differs between "Test Case
 *        1" (purely diffusive, no velocity-dependent dispersion) and Test Cases 2/3
 *        (velocity-dependent dispersion via spatialparams.hh's dispersionAlphas(),
 *        smaller molecular diffusion).
 */
#ifndef DUMUX_TEST_HENRY_FAHS_FLUIDSYSTEM_HH
#define DUMUX_TEST_HENRY_FAHS_FLUIDSYSTEM_HH

#include <string>
#include <limits>

#include <dumux/io/name.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidsystems/base.hh>

namespace Dumux::FluidSystems {

/*!
 * \ingroup OnePNCTests
 * \brief Single-phase, two-component fluid system matching Fahs et al. (2016), Table
 *        2: rho0 = 1000 kg/m^3, rho1 = 1025 kg/m^3, mu = 1e-3 Pa s (identical across
 *        all three test cases); the molecular diffusion coefficient Dm is read from
 *        Problem.MolecularDiffusion (18.86e-6 m^2/s for Test Case 1, 9.43e-8 m^2/s for
 *        Test Cases 2 and 3).
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
    //! Molecular diffusion coefficient [m^2/s], read from Problem.MolecularDiffusion
    //! (see the class documentation above for the value used in each test case)
    template<class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx, int compIIdx, int compJIdx)
    {
        static const Scalar Dm = getParam<Scalar>("Problem.MolecularDiffusion");
        return Dm;
    }
};

} // end namespace Dumux::FluidSystems

#endif
