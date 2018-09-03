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
 * \ingroup ThreePWaterOilModel
 * \brief Adds I/O fields specific to the twop model
 */
#ifndef DUMUX_3P2CNI_IO_FIELDS_HH
#define DUMUX_3P2CNI_IO_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilModel
 * \brief Adds I/O fields specific to the three-phase three-component model
 */
template <class Indices>
class ThreePWaterOilIOFields
{

public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized out output fields
        out.addVolumeVariable( [](const auto& v){ return v.saturation(FluidSystem::wPhaseIdx); }, "S_w");
        out.addVolumeVariable( [](const auto& v){ return v.saturation(FluidSystem::nPhaseIdx); },"S_n");
        out.addVolumeVariable( [](const auto& v){ return v.saturation(FluidSystem::gPhaseIdx); },"S_g");
        out.addVolumeVariable( [](const auto& v){ return v.pressure(FluidSystem::wPhaseIdx); },"p_w");
        out.addVolumeVariable( [](const auto& v){ return v.pressure(FluidSystem::nPhaseIdx); },"p_n");
        out.addVolumeVariable( [](const auto& v){ return v.pressure(FluidSystem::gPhaseIdx); },"p_g");
        out.addVolumeVariable( [](const auto& v){ return v.density(FluidSystem::wPhaseIdx); },"rho_w");
        out.addVolumeVariable( [](const auto& v){ return v.density(FluidSystem::nPhaseIdx); },"rho_n");
        out.addVolumeVariable( [](const auto& v){ return v.density(FluidSystem::gPhaseIdx); },"rho_g");
        out.addVolumeVariable( [](const auto& v){ return v.mobility(FluidSystem::wPhaseIdx); },"mob_w");
        out.addVolumeVariable( [](const auto& v){ return v.mobility(FluidSystem::nPhaseIdx); },"mob_n");
        out.addVolumeVariable( [](const auto& v){ return v.mobility(FluidSystem::gPhaseIdx); },"mob_g");
        out.addVolumeVariable( [](const auto& v){ return v.viscosity(FluidSystem::wPhaseIdx); },"viscos_w");
        out.addVolumeVariable( [](const auto& v){ return v.viscosity(FluidSystem::nPhaseIdx); },"viscos_n");
        out.addVolumeVariable( [](const auto& v){ return v.viscosity(FluidSystem::gPhaseIdx); },"viscos_g");

        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
            for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                out.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },
                                      "x^" + FluidSystem::componentName(j) + "_" +  FluidSystem::phaseName(i));

        out.addVolumeVariable( [](const auto& v){ return v.porosity(); },"porosity");
        out.addVolumeVariable( [](const auto& v){ return v.priVars().state(); },"phase presence");
        out.addVolumeVariable( [](const auto& v){ return v.permeability(); },"permeability");
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
        switch (state)
        {
            case Indices::threePhases:
            {
                const std::vector<std::string> s1 = {"p_g",
                                                     "S_w",
                                                     "S_n"};
                return s1[pvIdx];
            }
            case Indices::wPhaseOnly:
            {
                const std::vector<std::string> s2 = {"p_w",
                                                     "T",
                                                     "x^" + FluidSystem::componentName(FluidSystem::nCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::wPhaseIdx)};
                return s2[pvIdx];
            }
            case Indices::gnPhaseOnly:
            {
                const std::vector<std::string> s3 = {"p_g",
                                                     "S_n",
                                                     "x^" + FluidSystem::componentName(FluidSystem::wCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::nPhaseIdx)};
                return s3[pvIdx];
            }
            case Indices::wnPhaseOnly:
            {
                const std::vector<std::string> s4 = {"p_w",
                                                     "T",
                                                     "S_n"};
                return s4[pvIdx];
            }
            case Indices::gPhaseOnly:
            {
                const std::vector<std::string> s5 = {"p_g",
                                                     "T",
                                                     "x^" + FluidSystem::componentName(FluidSystem::nCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::gPhaseIdx)};
                return s5[pvIdx];
            }
            case Indices::wgPhaseOnly:
            {
                const std::vector<std::string> s6 = {"p_g",
                                                     "S_w",
                                                     "x^" + FluidSystem::componentName(FluidSystem::nCompIdx) + "_" +  FluidSystem::phaseName(FluidSystem::gPhaseIdx)};
                return s6[pvIdx];
            }
        }
    }
};

} // end namespace Dumux

#endif
