// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup InputOutput
 * \brief A collection of input/output field names for common physical quantities
 */
#ifndef DUMUX_IO_FIELD_NAMES_HH
#define DUMUX_IO_FIELD_NAMES_HH

#include <string>

namespace Dumux {
namespace IOFieldNames {

//! name of variable pressure
template<class FluidSystem>
std::string pressure(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "p" : "p_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable saturation
template<class FluidSystem>
std::string saturation(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "S" : "S_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable temperature (equilibrium models)
std::string temperature() noexcept
{ return "T"; }

//! name of variable temperature (non-equilibrium models)
template<class FluidSystem>
std::string fluidTemperature(int phaseIdx) noexcept
{ return "T_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable temperature (non-equilibrium models)
std::string solidTemperature() noexcept
{ return "T_s"; }

//! name of variable density
template<class FluidSystem>
std::string density(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "rho" : "rho_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable molar density
template<class FluidSystem>
std::string molarDensity(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "rhoMolar" : "rhoMolar_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable mobility
template<class FluidSystem>
std::string mobility(int phaseIdx) noexcept
{ return (FluidSystem::numPhases == 1) ? "mob" : "mob_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable mole fraction
template<class FluidSystem>
std::string moleFraction(int phaseIdx, int compIdx) noexcept
{ return "x^" + FluidSystem::componentName(compIdx) + "_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable mass fraction
template<class FluidSystem>
std::string moleFraction(int phaseIdx, int compIdx) noexcept
{ return "X^" + FluidSystem::componentName(compIdx) + "_" + FluidSystem::phaseName(phaseIdx); }

//! name of variable capillary pressure
std::string capillaryPressure() noexcept
{ return "pc"; }

//! name of variable porosity
std::string porosity() noexcept
{ return "porosity"; }

//! name of variable phase presence
std::string phasePresence() noexcept
{ return "phase presence"; }

//! name of solid volume fraction
template<class SolidSystem>
std::string solidVolumeFraction(int compIdx) noexcept
{ return "precipitateVolumeFraction^" + SolidSystem::componentName(compIdx); }

} // end namespace IO
} // end namespace Dumux

#endif
