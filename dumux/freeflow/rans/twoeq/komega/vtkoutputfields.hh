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
 * \ingroup KOmegaModel
 * \copydoc Dumux::KOmegaVtkOutputFields
 */
#ifndef DUMUX_KOMEGA_VTK_OUTPUT_FIELDS_HH
#define DUMUX_KOMEGA_VTK_OUTPUT_FIELDS_HH

#include <dumux/freeflow/rans/vtkoutputfields.hh>

namespace Dumux
{

/*!
 * \ingroup KOmegaModel
 * \brief Adds vtk output fields for the Reynolds-Averaged Navier-Stokes model
 */
template<class FVGridGeometry>
class KOmegaVtkOutputFields : public RANSVtkOutputFields<FVGridGeometry>
{
    enum { dim = FVGridGeometry::GridView::dimension };

public:
    //! Initialize the Navier-Stokes specific vtk output fields.
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        RANSVtkOutputFields<FVGridGeometry>::init(vtk);
        add(vtk);
    }

    //! Add the KOmegaModel specific vtk output fields.
    template <class VtkOutputModule>
    static void add(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const auto& v){ return v.turbulentKineticEnergy(); }, "turbulentKineticEnergy");
        vtk.addVolumeVariable([](const auto& v){ return v.dissipation(); }, "dissipation");
        vtk.addVolumeVariable([](const auto& v){ return v.stressTensorScalarProduct(); }, "stressTensorScalarProduct");
        vtk.addVolumeVariable([](const auto& v){ return v.kinematicEddyViscosity(); }, "eddyViscosity");
        vtk.addVolumeVariable([](const auto& v){ return 2.0 * v.kinematicEddyViscosity() * v.stressTensorScalarProduct();}, "production_k");
        vtk.addVolumeVariable([](const auto& v){ return 2.0 * v.kinematicEddyViscosity() * v.stressTensorScalarProduct() * v.alpha() * ( v.dissipation() / v.turbulentKineticEnergy() ) ; }, "production_omega");
        vtk.addVolumeVariable([](const auto& v){ return v.betaK() * v.turbulentKineticEnergy() * v.dissipation() ;}, "destruction_k");
        vtk.addVolumeVariable([](const auto& v){ return v.betaOmega() * v.dissipation() * v.dissipation() ;}, "destruction_omega");
    }
};

} // end namespace Dumux

#endif
