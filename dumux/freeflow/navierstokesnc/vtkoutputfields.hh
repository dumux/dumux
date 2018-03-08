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
 * \ingroup NavierStokesNCModel
 * \copydoc Dumux::NavierStokesNCVtkOutputFields
 */
#ifndef DUMUX_NAVIER_STOKES_NC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_NAVIER_STOKES_NC_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/vtkoutputfields.hh>

namespace Dumux
{

/*!
 * \ingroup NavierStokesNCModel
 * \brief Adds vtk output fields specific to the NavierStokesNC model
 */
template<class FVGridGeometry, class FluidSystem, int phaseIdx>
class NavierStokesNCVtkOutputFields
{

public:
    //! Initialize the Navier-StokesNC specific vtk output fields.
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        NavierStokesVtkOutputFields<FVGridGeometry>::init(vtk);

        vtk.addVolumeVariable([](const auto& v){ return v.molarDensity(); }, "rhoMolar");
        vtk.addVolumeVariable([](const auto& v){ return v.density(); }, "rho");

        for (int j = 0; j < FluidSystem::numComponents; ++j)
        {
            vtk.addVolumeVariable([j](const auto& v){ return v.massFraction(phaseIdx,j); }, "X^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(phaseIdx));
            vtk.addVolumeVariable([j](const auto& v){ return v.moleFraction(phaseIdx,j); }, "x^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(phaseIdx));
        }
    }
};

} // end namespace Dumux

#endif
