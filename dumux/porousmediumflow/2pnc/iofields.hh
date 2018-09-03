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
 * \ingroup TwoPNCModel
 * \brief Adds I/O fields specific to the twop-nc model
 */
#ifndef DUMUX_TWOP_NC_IO_FIELDS_HH
#define DUMUX_TWOP_NC_IO_FIELDS_HH

#include <dumux/porousmediumflow/2p/iofields.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNCModel
 * \brief Adds I/O fields specific to the TwoPNC model
 */
template <TwoPFormulation priVarFormulation, class Indices, bool useMoles, bool setMoleFractionsForFirstPhase>
class TwoPNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // use default fields from the 2p model
        TwoPIOFields<priVarFormulation>::initOutputModule(out);

        //output additional to TwoP output:
        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
            for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                out.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },
                                    "x^"+ FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(i));

        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
            out.addVolumeVariable([i](const auto& v){ return v.molarDensity(i); },
                                    "rhoMolar_" + FluidSystem::phaseName(i));

        if (VolumeVariables::numComponents() < 3){
            for (int i = 0; i < VolumeVariables::numPhases(); ++i)
                for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                    out.addVolumeVariable([i,j](const auto& v){ return v.massFraction(i,j); },
                                    "X^"+ FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(i));
        }

        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); }, "phase presence");
    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        const std::string xString = useMoles ? "x" : "X";

        std::string phaseNameSecComps;
        if (state == Indices::firstPhaseOnly
            || (state == Indices::bothPhases && setMoleFractionsForFirstPhase))
            phaseNameSecComps = FluidSystem::phaseName(FluidSystem::phase0Idx);
        else
            phaseNameSecComps = FluidSystem::phaseName(FluidSystem::phase1Idx);

        if (pvIdx > 1)
            return xString + "^" + FluidSystem::componentName(pvIdx) + "_" + phaseNameSecComps;

        const std::vector<std::string> p0s1SwitchedPvNames = {
            xString + "^" + FluidSystem::componentName(FluidSystem::comp1Idx) + "_" + FluidSystem::phaseName(FluidSystem::phase0Idx),
            xString + "^" + FluidSystem::componentName(FluidSystem::comp0Idx) + "_" + FluidSystem::phaseName(FluidSystem::phase1Idx),
            "S_n"};
        const std::vector<std::string> p1s0SwitchedPvNames = {
            xString + "^" + FluidSystem::componentName(FluidSystem::comp1Idx) + "_" + FluidSystem::phaseName(FluidSystem::phase0Idx),
            xString + "^" + FluidSystem::componentName(FluidSystem::comp0Idx) + "_" + FluidSystem::phaseName(FluidSystem::phase1Idx),
            "S_w"};

        switch (priVarFormulation)
        {
        case TwoPFormulation::p0s1:
            return pvIdx == 0 ? "p_w" : p0s1SwitchedPvNames[state-1];
        case TwoPFormulation::p1s0:
            return pvIdx == 0 ? "p_n" : p1s0SwitchedPvNames[state-1];
        default: DUNE_THROW(Dune::InvalidStateException, "Invalid formulation ");
        }
    }
};


} // end namespace Dumux

#endif
