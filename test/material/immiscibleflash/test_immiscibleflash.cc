// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief This is a program to test the flash calculation which uses
 *        non-linear complementarity problems (NCP).
 *
 * A flash calculation determines the pressures, saturations and
 * composition of all phases given the total mass (or, as in this case
 * the total number of moles) in a given amount of pore space.
 */

#include <config.h>

#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/immiscibleflash.hh>

#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/material/fluidsystems/h2on2.hh>

#include <dumux/material/fluidmatrixinteractions/mp/mpadapter.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

#include <dune/common/exceptions.hh>

template <class Scalar, class FluidState>
void checkSame(const FluidState &fsRef, const FluidState &fsFlash)
{
    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        Scalar error;

        // check the pressures
        using std::abs;
        error = 1 - fsRef.pressure(phaseIdx)/fsFlash.pressure(phaseIdx);
        if (abs(error) > 1e-6) {
            DUNE_THROW(Dune::Exception, "Mismatch to reference. Pressure of phase "
            << phaseIdx << " calculated by flash is: "
            << fsFlash.pressure(phaseIdx)  << " vs "
            << fsRef.pressure(phaseIdx) << " calculated as reference."
            << " Error = " << error << "\n");
        }

        // check the saturations
        error = fsRef.saturation(phaseIdx) - fsFlash.saturation(phaseIdx);
        if (abs(error) > 1e-6) {
            DUNE_THROW(Dune::Exception, "Mismatch to reference. Saturation of phase "
            << phaseIdx << " calculated by flash is: "
            << fsFlash.saturation(phaseIdx)  << " vs "
            << fsRef.saturation(phaseIdx) << " calculated as reference."
            << " Error = " << error << "\n");
        }


        // check the compositions
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            error = fsRef.moleFraction(phaseIdx, compIdx) - fsFlash.moleFraction(phaseIdx, compIdx);
            if (abs(error) > 1e-6) {
                DUNE_THROW(Dune::Exception, "Mismatch to reference. Mole fraction of component "
                << compIdx << " in phase "
                << phaseIdx << " calculated by flash is: "
                << fsFlash.moleFraction(phaseIdx, compIdx) << " vs "
                << fsRef.moleFraction(phaseIdx, compIdx) << " calculated as reference."
                << " Error = " << error << "\n");
            }
        }
    }
}

template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void checkImmiscibleFlash(const FluidState &fsRef,
                          const MaterialLaw& material)
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

    // calculate the total amount of stuff in the reference fluid
    // phase
    ComponentVector globalMolarities(0.0);
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            globalMolarities[compIdx] +=
                fsRef.saturation(phaseIdx)*fsRef.molarity(phaseIdx, compIdx);
        }
    }

    // create the flash
    Dumux::ImmiscibleFlash<Scalar, FluidSystem> flash(/*wettingPhaseIdx=*/0);
    // initialize the fluid state for the flash calculation
    FluidState fsFlash;

    fsFlash.setTemperature(fsRef.temperature(/*phaseIdx=*/0));

    // run the flash calculation
    typename FluidSystem::ParameterCache paramCache;
    flash.guessInitial(fsFlash, paramCache, globalMolarities);
    flash.solve(fsFlash, paramCache, material, globalMolarities);

    // compare the "flashed" fluid state with the reference one
    checkSame<Scalar>(fsRef, fsFlash);
}


template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void completeReferenceFluidState(FluidState &fs,
                                 const MaterialLaw& material,
                                 int refPhaseIdx)
{
    enum { numPhases = FluidSystem::numPhases };

    int otherPhaseIdx = 1 - refPhaseIdx;

    // calculate the other saturation
    fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

    // calculate the capillary pressure
    const auto pc = material.capillaryPressures(fs, /*wPhaseIdx=*/0);
    fs.setPressure(otherPhaseIdx,
                   fs.pressure(refPhaseIdx)
                   + (pc[otherPhaseIdx] - pc[refPhaseIdx]));

    // set all phase densities
    typename FluidSystem::ParameterCache paramCache;
    paramCache.updateAll(fs);
    for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
        Scalar rhoMolar = FluidSystem::molarDensity(fs, paramCache, phaseIdx);
        fs.setDensity(phaseIdx, rho);
        fs.setMolarDensity(phaseIdx, rhoMolar);
    }
}


int main()
{
    using Scalar = double;
    using FluidSystem = Dumux::FluidSystems::H2ON2<Scalar>;
    using ImmiscibleFluidState = Dumux::ImmiscibleFluidState<Scalar, FluidSystem>;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    using PcKrSw = Dumux::FluidMatrix::BrooksCoreyDefault<Scalar>;

    Scalar T = 273.15 + 25;

    // initialize the tables of the fluid system
    Scalar Tmin = T - 1.0;
    Scalar Tmax = T + 1.0;
    int nT = 3;

    Scalar pmin = 0.0;
    Scalar pmax = 1.25 * 2e6;
    int np = 100;

    FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

    // set the parameters for the capillary pressure law
    typename PcKrSw::BasicParams bParams(/*pe*/0.0, /*lambda*/2.0);
    auto pcKrSw = PcKrSw(bParams);
    auto material = Dumux::FluidMatrix::MPAdapter(pcKrSw);

    ImmiscibleFluidState fsRef;

    // create an fluid state which is consistent

    // set the fluid temperatures
    fsRef.setTemperature(T);

    ////////////////
    // only liquid
    ////////////////
    std::cout << "testing single-phase liquid\n";

    // set liquid saturation and pressure
    fsRef.setSaturation(liquidPhaseIdx, 1.0);
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem>(fsRef, material, liquidPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem>(fsRef, material);

    ////////////////
    // only gas
    ////////////////
    std::cout << "testing single-phase gas\n";

    // set gas saturation and pressure
    fsRef.setSaturation(gasPhaseIdx, 1.0);
    fsRef.setPressure(gasPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem>(fsRef, material, gasPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem>(fsRef, material);

    ////////////////
    // both phases
    ////////////////
    std::cout << "testing two-phase\n";

    // set liquid saturation and pressure
    fsRef.setSaturation(liquidPhaseIdx, 0.5);
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem>(fsRef, material, liquidPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem>(fsRef, material);

    ////////////////
    // with capillary pressure
    ////////////////
    std::cout << "testing two-phase with capillary pressure\n";

    typename PcKrSw::BasicParams bParams2(/*pe*/1e3, /*lambda*/2.0);
    auto pcKrSw2 = PcKrSw(bParams2);
    auto material2 = Dumux::FluidMatrix::MPAdapter(pcKrSw2);

    // set liquid saturation
    fsRef.setSaturation(liquidPhaseIdx, 0.5);

    // set pressure of the liquid phase
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // set the remaining parameters of the reference fluid state
    completeReferenceFluidState<Scalar, FluidSystem>(fsRef, material2, liquidPhaseIdx);

    // check the flash calculation
    checkImmiscibleFlash<Scalar, FluidSystem>(fsRef, material2);

    return 0;
}
