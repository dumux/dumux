// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief A collection of input/output field names for common physical quantities
 */
#ifndef DUMUX_IO_NAME_HH
#define DUMUX_IO_NAME_HH

#include <string>

namespace Dumux {
namespace IOName {

//! I/O name of pressure for multiphase systems
template<class FluidSystem>
std::string pressure(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "p" : "p_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of pressure for singlephase systems
std::string pressure() noexcept
{ return "p"; }

//! I/O name of saturation for multiphase systems
template<class FluidSystem>
std::string saturation(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "S" : "S_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of saturation for singlephase systems
std::string saturation() noexcept
{ return "S"; }

//! I/O name of temperature for equilibrium models
std::string temperature() noexcept
{ return "T"; }

//! I/O name of temperature for non-equilibrium models
template<class FluidSystem>
std::string fluidTemperature(int phaseIdx) noexcept
{ return "T_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of solid temperature for non-equilibrium models
std::string solidTemperature() noexcept
{ return "T_s"; }

//! I/O name of density for multiphase systems
template<class FluidSystem>
std::string density(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "rho" : "rho_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of density for singlephase systems
std::string density() noexcept
{ return "rho"; }

//! I/O name of viscosity for multiphase systems
template<class FluidSystem>
std::string viscosity(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "mu" : "mu_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of viscosity for singlephase systems
std::string viscosity() noexcept
{ return "mu"; }

//! I/O name of molar density for multiphase systems
template<class FluidSystem>
std::string molarDensity(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "rhoMolar" : "rhoMolar_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of molar density for singlephase systems
std::string molarDensity() noexcept
{ return "rhoMolar"; }

//! I/O name of relative permeability for multiphase systems
template<class FluidSystem>
std::string relativePermeability(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "kr" : "kr_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of relative permeability for singlephase systems
std::string relativePermeability() noexcept
{ return "kr"; }

//! I/O name of mobility for multiphase systems
template<class FluidSystem>
std::string mobility(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "mob" : "mob_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of mobility for singlephase systems
std::string mobility() noexcept
{ return "mob"; }

//! I/O name of mole fraction
template<class FluidSystem>
std::string moleFraction(int phaseIdx, int compIdx) noexcept
{ return "x^" + FluidSystem::componentName(compIdx) + "_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of mass fraction
template<class FluidSystem>
std::string massFraction(int phaseIdx, int compIdx) noexcept
{ return "X^" + FluidSystem::componentName(compIdx) + "_" + FluidSystem::phaseName(phaseIdx); }

//! I/O name of liquid phase
std::string liquidPhase() noexcept
{ return "liq"; }

//! I/O name of gaseous phase
std::string gaseousPhase() noexcept
{ return "gas"; }

//! I/O name of aqueous phase
std::string aqueousPhase() noexcept
{ return "aq"; }

//! I/O name of napl phase
std::string naplPhase() noexcept
{ return "napl"; }

//! I/O name of capillary pressure
std::string capillaryPressure() noexcept
{ return "pc"; }

//! I/O name of porosity
std::string porosity() noexcept
{ return "porosity"; }

//! I/O name of permeability
std::string permeability() noexcept
{ return "permeability"; }

//! I/O name of phase presence
std::string phasePresence() noexcept
{ return "phase presence"; }

//! I/O name of pressure head
std::string pressureHead() noexcept
{ return "pressure head"; }

//! I/O name of water content
std::string waterContent() noexcept
{ return "water content"; }

//! I/O name of solid volume fraction
template<class SolidSystem>
std::string solidVolumeFraction(int compIdx = 0) noexcept
{ return "precipitateVolumeFraction^" + SolidSystem::componentName(compIdx); }

//! I/O name of displacement
std::string displacement() noexcept
{ return "u"; }

} // end namespace IOName
} // end namespace Dumux

#endif
