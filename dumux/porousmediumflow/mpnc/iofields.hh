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
 * \ingroup MPNCModel
 * \brief Adds I/O fields specific to the mpnc model.
 */

#ifndef DUMUX_MPNC_IO_FIELDS_HH
#define DUMUX_MPNC_IO_FIELDS_HH

#include <dune/common/exceptions.hh>

#include <dumux/io/name.hh>
#include "pressureformulation.hh"

namespace Dumux {

/*!
 * \ingroup MPNCModel
 * \brief Adds I/O fields specific to the mpnc model.
 */
class MPNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        for (int i = 0; i < VolumeVariables::numFluidPhases(); ++i)
        {
            out.addVolumeVariable([i](const auto& v){ return v.saturation(i); },
                                  IOName::saturation<FluidSystem>(i));
            out.addVolumeVariable([i](const auto& v){ return v.pressure(i); },
                                  IOName::pressure<FluidSystem>(i));
            out.addVolumeVariable([i](const auto& v){ return v.density(i); },
                                  IOName::density<FluidSystem>(i));
            out.addVolumeVariable([i](const auto& v){ return v.mobility(i); },
                                  IOName::mobility<FluidSystem>(i));

            for (int j = 0; j < VolumeVariables::numFluidComponents(); ++j)
                out.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },
                                      IOName::moleFraction<FluidSystem>(i, j));
        }

        for (int j = 0; j < VolumeVariables::numFluidComponents(); ++j)
        out.addVolumeVariable([j](const auto& v){ return v.fugacity(j); },
                              "fugacity^"+ FluidSystem::componentName(j));


        out.addVolumeVariable([](const auto& v){ return v.porosity(); },
                              IOName::porosity());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state=0)
    {
        if (pvIdx < ModelTraits::numFluidComponents())
            return "fugacity^"+ FluidSystem::componentName(pvIdx);
        else if (pvIdx < ModelTraits::numEq()-1)
            return "S_"+ FluidSystem::phaseName(pvIdx - ModelTraits::numFluidComponents());
        else
        {
            switch (ModelTraits::pressureFormulation())
            {
                case MpNcPressureFormulation::mostWettingFirst :
                    return "p_"+ FluidSystem::phaseName(0);
                case MpNcPressureFormulation::leastWettingFirst :
                    return "p_"+ FluidSystem::phaseName(ModelTraits::numFluidPhases()-1);
                default: DUNE_THROW(Dune::InvalidStateException, "Invalid formulation ");
            }
        }
    }

};

} // end namespace Dumux

#endif
