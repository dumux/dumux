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
 * \ingroup FreeflowNCModel
 * \copydoc Dumux::FreeflowNCVtkOutputFields
 */
#ifndef DUMUX_FREEFLOW_NC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_FREEFLOW_NC_VTK_OUTPUT_FIELDS_HH

#include <dumux/freeflow/navierstokes/vtkoutputfields.hh>

namespace Dumux
{

/*!
 * \ingroup FreeflowNCModel
 * \brief Adds vtk output fields specific to the FreeflowNC model
 */
template<class BaseVtkOutputFields, class ModelTraits, class FVGridGeometry, class FluidSystem>
class FreeflowNCVtkOutputFields
{

public:
    //! Initialize the FreeflowNC specific vtk output fields.
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        BaseVtkOutputFields::init(vtk);
        add(vtk);
    }

    //! Add the FreeflowNC specific vtk output fields.
    template <class VtkOutputModule>
    static void add(VtkOutputModule& vtk)
    {
        for (int j = 0; j < FluidSystem::numComponents; ++j)
        {
            vtk.addVolumeVariable([j](const auto& v){ return v.massFraction(j); }, "X^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(0));
            vtk.addVolumeVariable([j](const auto& v){ return v.moleFraction(j); }, "x^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(0));

            if (j != FluidSystem::getMainComponent(0))
            {
                vtk.addVolumeVariable([j](const auto& v){ return v.diffusionCoefficient(0, j); }, "D^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(0));

                // the eddy diffusivity is recalculated for an arbitrary component which is not the phase component
                if (ModelTraits::usesTurbulenceModel())
                    vtk.addVolumeVariable([j](const auto& v){ return v.effectiveDiffusivity(0, j) - v.diffusionCoefficient(0, j); }, "D_t");
            }
        }
    }
};

} // end namespace Dumux

#endif
