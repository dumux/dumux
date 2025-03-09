// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FrictionLaws
 * \brief A pseudo friction law with no bottom friction
 */

#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NOFRICTION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NOFRICTION_HH

#include "frictionlaw.hh"

namespace Dumux {
/*!
 * \addtogroup FrictionLaws
 * \copydetails Dumux::FrictionLawNoFriction
 */

/*!
 * \ingroup FrictionLaws
 * \brief A pseudo friction law with no bottom friction
 *
 * ### No Friction
 *
 * This friction law sets the stress between the flowing fluid and the bottom,
 * which is called bottom shear stress, to zero.
 * The bottom shear stress is needed to calculate on the one hand the loss of
 * momentum due to bottom friction and on the other hand the bedload transport rate.
 *
 */
template <typename VolumeVariables>
class FrictionLawNoFriction : public FrictionLaw<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     */
    FrictionLawNoFriction() = default;

    /*!
     * \brief Compute the bottom shear stress.
     *
     * \param volVars Volume variables
     *
     * Compute the bottom shear stress due to bottom friction.
     * The bottom shear stress is a projection of the shear stress tensor onto the river bed.
     * It can therefore be represented by a (tangent) vector with two entries.
     * For this law without bottom friction, the bottom shear stress is zero.
     *
     * \return shear stress in N/m^2. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> bottomShearStress(const VolumeVariables& volVars) const final
    {
        return {0.0, 0.0};
    }
};

} // end namespace Dumux

#endif
