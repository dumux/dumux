// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowNCModel
 * \copydoc Dumux::NavierStokesMassOnePNCIOFields
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1PNC_IO_FIELDS_HH
#define DUMUX_NAVIERSTOKES_MASS_1PNC_IO_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief Adds I/O fields specific to the FreeflowNC model
 */
template<class BaseOutputFields>
struct NavierStokesMassOnePNCIOFields
{
    //! Initialize the FreeflowNC specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        BaseOutputFields::initOutputModule(out);

        // TODO only output molar density if we use mole fractions
        out.addVolumeVariable([](const auto& v){ return v.molarDensity(); }, IOName::molarDensity());

        using FluidSystem = typename OutputModule::VolumeVariables::FluidSystem;
        for (int j = 0; j < FluidSystem::numComponents; ++j)
        {
            out.addVolumeVariable([j](const auto& v){ return v.massFraction(j); }, IOName::massFraction<FluidSystem>(0, j));
            out.addVolumeVariable([j](const auto& v){ return v.moleFraction(j); }, IOName::moleFraction<FluidSystem>(0, j));

            if (j != FluidSystem::getMainComponent(0))
            {
                out.addVolumeVariable([j](const auto& v){ return v.diffusionCoefficient(0,0, j); }, "D^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(0));
            }
        }
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        // priVars: p, x_0, ..., x_numComp-1, otherPv ..., T
        if (pvIdx > 0 && pvIdx < ModelTraits::numFluidComponents() + 1)
            return ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(0, pvIdx - 1)
                                           : IOName::massFraction<FluidSystem>(0, pvIdx - 1);
        else
            return BaseOutputFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx, state);

    }

    template<class VolumeVariables>
    static double getEffectiveDiffusionCoefficient_(const VolumeVariables& volVars, const int phaseIdx, const int compIdx)
    {
        return volVars.effectiveDiffusionCoefficient(phaseIdx, VolumeVariables::FluidSystem::getMainComponent(phaseIdx), compIdx);
    }
};

} // end namespace Dumux

#endif
