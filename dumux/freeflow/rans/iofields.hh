// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

        // velocityGradients is a tensor, gradient of each velocity component is added here.
        out.addVolumeVariable([](const auto& v){ return v.velocityGradients()[0]; }, "dv_x/ds_");
        if (dim > 1)
            out.addVolumeVariable([](const auto& v){ return v.velocityGradients()[1]; }, "dv_y/ds_");
        if (dim > 2)
            out.addVolumeVariable([](const auto& v){ return v.velocityGradients()[2]; }, "dv_z/ds_");
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
