// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FrictionLaws
 * \brief Implementation of the friction law after Manning.
 */

#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_MANNING_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_MANNING_HH

#include <algorithm>
#include <cmath>
#include <dune/common/math.hh>
#include "frictionlaw.hh"

namespace Dumux {
/*!
 * \addtogroup FrictionLaws
 * \copydetails Dumux::FrictionLawManning
 */

/*!
 * \ingroup FrictionLaws
 * \brief Implementation of the friction law after Manning.
 *
 * ### Manning
 *
 * This friction law calculates the stress between the flowing fluid and the bottom,
 * which is called bottom shear stress, using the Manning friction law:
 *
 * \f$\tau_{x} = \frac{g}{(\frac{h^{1/6}}{n})^2} u \sqrt{u^2 + v^2}\f$ and
 * \f$\tau_{y} = \frac{g}{(\frac{h^{1/6}}{n})^2} v \sqrt{u^2 + v^2}\f$
 *
 * with the gravity constant \f$\mathrm{g}\f$ in \f$\mathrm{[m/s^2]}\f$, the water depth
 * \f$\mathrm{h}\f$ in \f$\mathrm{[m]}\f$ and the Manning friction coefficient
 * \f$\mathrm{n}\f$ in \f$\mathrm{[s/m^{1/3}]}\f$.
 *
 * The bottom shear stress is needed to calculate on the one hand the loss of
 * momentum due to bottom friction and on the other hand the bedload transport rate.
 *
 * The LET mobility model is used to limit the friction for small water
 * depths if a roughness height > 0.0 is provided (default roughnessHeight = 0.0).
 */

template <typename VolumeVariables>
class FrictionLawManning : public FrictionLaw<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     *
     * \param gravity Gravity constant (in m/s^2)
     * \param manningN Manning friction coefficient (in s/m^(1/3)
     * \param roughnessHeight roughness height for limiting (in m) default = 0.0
     */
    FrictionLawManning(const Scalar gravity, const Scalar manningN, const Scalar roughnessHeight=0.0)
    : gravity_(gravity), manningN_(manningN), roughnessHeight_(roughnessHeight) {}

    /*!
     * \brief Compute the bottom shear stress.
     *
     * \param volVars Volume variables
     *
     * Compute the bottom shear stress due to bottom friction.
     * The bottom shear stress is a projection of the shear stress tensor onto the river bed.
     * It can therefore be represented by a (tangent) vector with two entries.
     *
     * \return shear stress in N/m^2. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> bottomShearStress(const VolumeVariables& volVars) const final
    {
        using std::pow;
        using Dune::power;
        using std::hypot;

        Dune::FieldVector<Scalar, 2> shearStress(0.0);

        const Scalar artificialWaterDepth = this->limitRoughH(roughnessHeight_, volVars.waterDepth());
        // c has units of m^(1/2)/s so c^2 has units of m/s^2
        const Scalar c = pow(volVars.waterDepth() + artificialWaterDepth, 1.0/6.0) * 1.0/(manningN_);
        const Scalar uv = hypot(volVars.velocity(0), volVars.velocity(1));
        const Scalar dimensionlessFactor = gravity_/(c*c);

        shearStress[0] = dimensionlessFactor * volVars.velocity(0) * uv * volVars.density();
        shearStress[1] = dimensionlessFactor * volVars.velocity(1) * uv * volVars.density();

        return shearStress;
    }

private:
    Scalar gravity_;
    Scalar manningN_;
    Scalar roughnessHeight_;
};

} // end namespace Dumux

#endif
