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
 * \ingroup MineralizationModel
 * \brief Adds vtk output fields specific to the models considering
 *        mineralization processes.
 */
#ifndef DUMUX_MINERALIZATION_VTK_OUTPUT_FIELDS_HH
#define DUMUX_MINERALIZATION_VTK_OUTPUT_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup MineralizationModel
 * \brief Adds vtk output fields specific to a NCMin model
 */
template<class NonMineralizationVtkOutputFields>
class MineralizationVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using SolidSystem = typename VtkOutputModule::VolumeVariables::SolidSystem;

        // output of the model without mineralization
        NonMineralizationVtkOutputFields::init(vtk);

        // additional output
        for (int i = 0; i < SolidSystem::numComponents - SolidSystem::numInertComponents; ++i)
        {
            vtk.addVolumeVariable([i](const auto& v){ return v.solidVolumeFraction(i); },"precipitateVolumeFraction_"+ SolidSystem::componentName(i));
        }
    }
};

} // end namespace Dumux

#endif
