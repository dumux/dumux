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
 * \ingroup OnePNCModel
 * \brief Adds vtk output fields specific to the OnePNC model
 */
#ifndef DUMUX_ONEPNC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_ONEPNC_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief Adds vtk output fields specific to the OnePNC model
 */
template<class TypeTag>
class OnePNCVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.pressure(phaseIdx); }, "pressure");
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.density(phaseIdx); }, "rho");
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.viscosity(phaseIdx); }, "mu");
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.pressure(phaseIdx) - 1e5; }, "delp");

        for (int i = 0; i < numComponents; ++i)
           vtk.addVolumeVariable([i](const auto& volVars){ return volVars.moleFraction(phaseIdx, i); },
                                     "x_" + std::string(FluidSystem::componentName(i)));

        for (int i = 0; i < numComponents; ++i)
           vtk.addVolumeVariable([i](const auto& volVars){ return volVars.massFraction(phaseIdx,i); },
                                     "X_" + std::string(FluidSystem::componentName(i)));
    }
};

} // end namespace Dumux

#endif
