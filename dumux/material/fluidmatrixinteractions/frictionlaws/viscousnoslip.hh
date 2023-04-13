// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \copydoc Dumux::FrictionLawViscousNoSlip
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_VISCOUS_NOSLIP_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_VISCOUS_NOSLIP_HH

#include <algorithm>
#include <cmath>
#include <dune/common/math.hh>

#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of a viscous no-slip bottom friction law
 *
 * This assumes thin film flow with a parabolic velocity profile in depth
 * (for the depth-averaged shallow water equations). The velocity profile
 * and associated bottom shear stress can be derived from plane Poiseuille flow
 * with a free surface boundary condition on top and a no-slip boundary condition
 * on the bottom.
 */

template <typename VolumeVariables>
class FrictionLawViscousNoSlip : public FrictionLaw<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Compute the bottom shear stress.
     *
     * Compute the bottom shear stress due to bottom friction.
     * The bottom shear stress is a projection of the shear stress tensor onto the bottom plane.
     * It can therefore be represented by a (tangent) vector with two entries.
     *
     * \return shear stress in N/m^2. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> bottomShearStress(const VolumeVariables& volVars) const final
    {
        // assume a parabolic velocity profile with no-slip BC on the bottom
        // and zero stress condition on the free surface
        // note that the velocity corresponds to the height-averaged velocity
        Dune::FieldVector<Scalar, 2> shearStress(0.0);
        shearStress[0] = volVars.viscosity()*volVars.velocity(0) * 3.0 / volVars.waterDepth();
        shearStress[1] = volVars.viscosity()*volVars.velocity(1) * 3.0 / volVars.waterDepth();
        return shearStress;
    }
};

} // end namespace Dumux

#endif
