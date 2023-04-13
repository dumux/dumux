// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \copydoc Dumux::FrictionLaw
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_HH

#include <algorithm>
#include <cmath>

#include <dune/common/fvector.hh>

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the abstract base class for friction laws.
 *
 * A LET mobility model can be used to add an artificial water depth to
 * limit the friction for small water depths.
 */

template <typename VolumeVariables >
class FrictionLaw
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Compute the bottom shear stress.
     *
     * \param volVars Volume variables
     *
     * Compute the bottom shear stress due to bottom friction.
     * The bottom shear stress is a projection of the shear stress tensor onto the river bed.
     * It can therefore be represented by a (tangent) vector with two entries.
     *
     * \return shear stress [N/m^2]. First entry is the x-component, the second the y-component.
     */
    virtual Dune::FieldVector<Scalar, 2> bottomShearStress(const VolumeVariables& volVars) const = 0;

    /*!
     * \brief Limit the friction for small water depth.
     *
     * Compute a small artificial water depth that is added to the
     * actual water depth to avoid extreme friction values which can
     * occur for small water depths.
     *
     * The function is called with a roughness height, which can be
     * seen as roughness height of the surface (e.g. grain size). For a
     * zero roughnessHeight the artificial water depth will be zero.
     *
     * A water depth minUpperH  (minUpperH = 2 * roughnessHeight) is
     * defined for the limiting. Limiting is applied between both
     * depths.
     *
     * ------------------------- minUpperH -----------
     *
     *
     *
     * ------------------------roughnessHeight ---------------
     *    /\  /\   roughness                  /grain\
     * -------------------------------bottom ------------------
     * /////////////////////////////////////////////////
     *
     * For the limiting the LET model is used, which is usually applied in the
     * porous media flow to limit the permeability due to the saturation. It employs
     * the three empirical parameters L, E and T, which describe the limiting curve (mobility).
     *
     * auto mobility = (mobility_max * pow(sw,L))/(pow(sw,L) + E * pow(1.0-sw,T));
     *
     * For the limitation of the roughness height L = 0.0, T = 2.0 and E = 1.0 are chosen.
     * Therefore the calculation of the mobility is simplified significantly.
     *
     * \param roughnessHeight roughness height of the representative structure (e.g. largest grain size).
     * \param waterDepth water depth.
     * \param eps If the roughness height falls below this threshold, this function returns zero.
     */
    Scalar limitRoughH(const Scalar roughnessHeight, const Scalar waterDepth, const Scalar eps=1.0e-12) const
    {
        using std::clamp;

        // return zero if the roughness depth is near zero
        if (roughnessHeight < eps) return 0.0;

        // calculate the artificial water depth
        const Scalar mobilityMax = 1.0; //!< maximal mobility

        const Scalar minUpperH = roughnessHeight * 2.0;
        const Scalar sw = clamp(waterDepth * (1.0/minUpperH), 0.0, 1.0);
        const Scalar mobility = mobilityMax /(1.0 + (1.0-sw)*(1.0-sw));
        return roughnessHeight * (1.0 - mobility);
    }

    // virtual base class needs a virtual destructor
    virtual ~FrictionLaw() = default;
};

} // end namespace Dumux

#endif
