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
 * \ingroup RANSModel
 * \copydoc Dumux::RANSIOFields
 */
#ifndef DUMUX_RANS_IO_FIELDS_HH
#define DUMUX_RANS_IO_FIELDS_HH

#include <dumux/freeflow/navierstokes/iofields.hh>

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Adds I/O fields for the Reynolds-Averaged Navier-Stokes model
 */
struct RANSIOFields
{
    //! Initialize the RANS specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        NavierStokesIOFields::initOutputModule(out);

        static const bool isFlatWallBounded = getParamFromGroup<bool>(out.paramGroup(), "RANS.IsFlatWallBounded", false);
        static const bool writeFlatWallBoundedFields = getParamFromGroup<bool>(out.paramGroup(), "RANS.WriteFlatWallBoundedFields", isFlatWallBounded);

        static constexpr auto dim = decltype(std::declval<typename OutputModule::VolumeVariables>().ccVelocityVector())::dimension;

        out.addVolumeVariable([](const auto& v){ return v.ccVelocityVector()[0] / v.velocityMaximum()[0]; }, "v_x/v_x,max");
        out.addVolumeVariable([](const auto& v){ return v.velocityGradients()[0]; }, "dv_x/dx_");
        if (dim > 1)
            out.addVolumeVariable([](const auto& v){ return v.velocityGradients()[1]; }, "dv_y/dx_");
        if (dim > 2)
            out.addVolumeVariable([](const auto& v){ return v.velocityGradients()[2]; }, "dv_z/dx_");
        out.addVolumeVariable([](const auto& v){ return v.pressure() - 1e5; }, "p_rel");
        out.addVolumeVariable([](const auto& v){ return v.viscosity() / v.density(); }, "nu");
        out.addVolumeVariable([](const auto& v){ return v.kinematicEddyViscosity(); }, "nu_t");
        out.addVolumeVariable([](const auto& v){ return v.wallDistance(); }, "l_w");
        if (writeFlatWallBoundedFields)
        {
            out.addVolumeVariable([](const auto& v){ return v.yPlus(); }, "y^+");
            out.addVolumeVariable([](const auto& v){ return v.uPlus(); }, "u^+");
        }
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        return NavierStokesIOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx, state);
    }
};

} // end namespace Dumux

#endif
