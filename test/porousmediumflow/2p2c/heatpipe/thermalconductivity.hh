// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Effective thermal conductivity of the heat-pipe solution of Huang et al. (2015).
 */
#ifndef DUMUX_HEATPIPE_THERMAL_CONDUCTIVITY_HH
#define DUMUX_HEATPIPE_THERMAL_CONDUCTIVITY_HH

#include <algorithm>
#include <cmath>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Effective thermal conductivity
 *        \f$\lambda = \lambda_\mathrm{dry} + \sqrt{S_w}(\lambda_\mathrm{wet} - \lambda_\mathrm{dry})\f$
 *        with the constant dry and wet conductivities of Huang et al. (2015).
 *
 * This has the form of ThermalConductivitySomertonTwoP, but with fixed end points
 * instead of geometric means of the solid and fluid conductivities.
 */
template<class Scalar>
class HeatPipeThermalConductivity
{
public:
    static constexpr Scalar lambdaWet = 1.13; // W/(m K)
    static constexpr Scalar lambdaDry = 0.582; // W/(m K)

    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;
        return effectiveThermalConductivity(volVars.saturation(FluidSystem::liquidPhaseIdx));
    }

    static Scalar effectiveThermalConductivity(const Scalar sw)
    {
        using std::max;
        using std::sqrt;
        return lambdaDry + sqrt(max(sw, Scalar(0.0)))*(lambdaWet - lambdaDry);
    }
};

} // end namespace Dumux

#endif
