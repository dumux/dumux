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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIOFields
 */
#ifndef DUMUX_NAVIER_STOKES_IO_FIELDS_HH
#define DUMUX_NAVIER_STOKES_IO_FIELDS_HH

#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/io/name.hh>

namespace Dumux
{

/*!
 * \ingroup NavierStokesModel
 * \brief Adds I/O fields for the Navier-Stokes model
 */
template<class FVGridGeometry>
class NavierStokesIOFields
{
    // Helper type used for tag dispatching (to add discretization-specific fields).
    template<DiscretizationMethod discMethod>
    using discMethodTag = std::integral_constant<DiscretizationMethod, discMethod>;

public:
    //! Initialize the Navier-Stokes specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using FluidSystem = typename OutputModule::VolumeVariables::FluidSystem;
        out.addVolumeVariable([](const auto& v){ return v.pressure(); }, IOName::pressure());
        out.addVolumeVariable([](const auto& v){ return v.molarDensity(); }, IOName::molarDensity<FluidSystem>());
        out.addVolumeVariable([](const auto& v){ return v.density(); }, IOName::density());

        // add discretization-specific fields
        additionalOutput_(out, discMethodTag<FVGridGeometry::discMethod>{});
    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    //! return the names of the primary variables
    template <class FluidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        const std::array<std::string, 3> velocities = {"v_x", "v_y", "v_z"};

        if (pvIdx < FVGridGeometry::Grid::dimension)
            return velocities[pvIdx];
        else
            return IOName::pressure();
    }

private:

    //! Adds discretization-specific fields (nothing by default).
    template <class OutputModule, class AnyMethod>
    static void additionalOutput_(OutputModule& out, AnyMethod)
    { }

    //! Adds discretization-specific fields (velocity vectors on the faces for the staggered discretization).
    template <class OutputModule>
    static void additionalOutput_(OutputModule& out, discMethodTag<DiscretizationMethod::staggered>)
    {
        const bool writeFaceVars = getParamFromGroup<bool>(out.paramGroup(), "Vtk.WriteFaceData", false);
        if(writeFaceVars)
        {
            auto faceVelocityVector = [](const typename FVGridGeometry::SubControlVolumeFace& scvf, const auto& faceVars)
                                      {
                                          using Scalar = typename OutputModule::VolumeVariables::PrimaryVariables::value_type;
                                          using VelocityVector = Dune::FieldVector<Scalar, FVGridGeometry::GridView::dimensionworld>;

                                          VelocityVector velocity(0.0);
                                          velocity[scvf.directionIndex()] = faceVars.velocitySelf();
                                          return velocity;
                                      };

            out.addFaceVariable(faceVelocityVector, "faceVelocity");

            auto faceNormalVelocity = [](const auto& faceVars)
                                      {
                                          return faceVars.velocitySelf();
                                      };

            out.addFaceVariable(faceNormalVelocity, "v");
        }
    }
};

} // end namespace Dumux

#endif
