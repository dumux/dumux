// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FrictionLaws
 * \brief Implementation of the abstract base class for friction laws.
 */

#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NIKURADSE_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NIKURADSE_HH

/*!
 * \file
 * \ingroup FrictionLaws
 * \brief Implementation of the friction law after Nikuradse.
 */
#include <algorithm>
#include <cmath>
#include <dune/common/math.hh>
#include "frictionlaw.hh"

namespace Dumux {
/*!
 * \addtogroup FrictionLaws
 * \copydetails Dumux::FrictionLawNikuradse
 */

/*!
 * \ingroup FrictionLaws
 * \brief Implementation of the friction law after Nikuradse.
 *
 * ### Nikuradse
 *
 * This friction law calculates the stress between the flowing fluid and the bottom,
 * which is called bottom shear stress, using the Nikuradse \cite Nikuradse1950 friction law
 *
 *\f$\tau_{x} = \frac{\kappa^2}{(ln\frac{12h}{ks})^2} u \sqrt{u^2 + v^2}\f$ and
 *\f$\tau_{y} = \frac{\kappa^2}{(ln\frac{12h}{ks})^2} v \sqrt{u^2 + v^2}\f$
 *
 * with the dimensionless Karman's constant \f$\mathrm{\kappa}\f$, the quivalent sand roughness
 * \f$\mathrm{ks}\f$ in \f$\mathrm{[m]}\f$ and the water depth \f$\mathrm{h}\f$
 * in \f$\mathrm{[m]}\f$.
 *
 * The bottom shear stress is needed to calculate on the one hand the loss of
 * momentum due to bottom friction and on the other hand the bedload transport rate.
 *
 * The LET mobility model is used to limit the friction for small water
 * depths if a roughness height > 0.0 is provided (default roughnessHeight = 0.0).
 *
 */

template <typename VolumeVariables>
class FrictionLawNikuradse : public FrictionLaw<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     *
     * \param ks Equivalent sand roughness (in m)
     * \param roughnessHeight roughness height (in m) default = 0.0
     */
    FrictionLawNikuradse(const Scalar ks, const Scalar roughnessHeight=0.0)
    : ks_(ks), roughnessHeight_(roughnessHeight) {}

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
        using Dune::power;
        using std::log;
        using std::hypot;

        Dune::FieldVector<Scalar, 2> shearStress(0.0);

        const Scalar artificialWaterDepth = this->limitRoughH(roughnessHeight_, volVars.waterDepth());
        const Scalar karmanConstant = 0.41; // Karman's constant is dimensionless
        const Scalar dimensionlessFactor = power(karmanConstant, 2)/power(log((12*(volVars.waterDepth() + artificialWaterDepth))/ks_), 2);
        const Scalar uv = hypot(volVars.velocity(0),volVars.velocity(1));

        shearStress[0] = dimensionlessFactor * volVars.velocity(0) * uv * volVars.density();
        shearStress[1] = dimensionlessFactor * volVars.velocity(1) * uv * volVars.density();

        return shearStress;
    }

private:
    Scalar ks_;
    Scalar roughnessHeight_;
};

} // end namespace Dumux

#endif
