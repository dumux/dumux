// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup FreeflowNCModel
 * \copydoc Dumux::FreeflowNCIOFields
 */
#ifndef DUMUX_FREEFLOW_NC_IO_FIELDS_HH
#define DUMUX_FREEFLOW_NC_IO_FIELDS_HH

#include <dumux/io/name.hh>
#include <dumux/freeflow/navierstokes/iofields.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief Adds I/O fields specific to the FreeflowNC model
 */
template<class BaseOutputFields, bool turbulenceModel = false>
struct FreeflowNCIOFields
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

                // the eddy diffusivity is recalculated for an arbitrary component which is not the phase component
                if (turbulenceModel)
                    out.addVolumeVariable([j](const auto& v){ return getEffectiveDiffusionCoefficient_(v, 0, j) - v.diffusionCoefficient(0,0, j); }, "D_t^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(0));
            }
        }
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        // priVars: v_0, ..., v_dim-1, p, x_0, ..., x_numComp-1, otherPv ..., T
        if (pvIdx > ModelTraits::dim() && pvIdx < ModelTraits::dim() + ModelTraits::numFluidComponents())
            return ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(0, pvIdx - ModelTraits::dim())
                                           : IOName::massFraction<FluidSystem>(0, pvIdx - ModelTraits::dim());
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
