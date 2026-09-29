// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \brief Regularization of the long-wave flux coefficient for vanishing depth and level water
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_REGULARIZATION_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_REGULARIZATION_HH

#include <algorithm>
#include <cmath>
#include <limits>

#include <dumux/common/parameters.hh>

namespace Dumux::LongWave {

/*!
 * \ingroup LongWaveModel
 * \brief The depth entering the hydraulic radius, bounded below by half the roughness height
 *
 * With the hydraulic radius equal to the depth, the conveyance \f$ h^{5/3} \f$ of a thin film has
 * a vanishing derivative at \f$ h = 0 \f$, so a dry degree of freedom has no entry on the diagonal
 * of the Jacobian. Below twice the roughness height `LongWave.Regularization.RoughnessHeight`
 * (default 1 mm) the hydraulic radius is smoothly (C1) bounded away from zero, which makes the
 * conveyance vanish linearly in the depth instead. A roughness height of zero disables this.
 */
template<class Scalar>
Scalar regularizedDepth(const Scalar h)
{
    using std::max, std::clamp;
    static const Scalar roughnessHeight = getParam<Scalar>("LongWave.Regularization.RoughnessHeight", 1e-3);
    if (roughnessHeight < 1e-20)
        return max(0.0, h);

    const Scalar minUpperH = roughnessHeight*2.0;
    const Scalar sw = clamp(h*(1.0/minUpperH), 0.0, 1.0);
    const Scalar mobility = 1.0/(1.0 + (1.0-sw)*(1.0-sw));
    return max(0.0, h) + roughnessHeight*(1.0 - mobility);
}

/*!
 * \ingroup LongWaveModel
 * \brief Peak of `regularInvSqrt` relative to \f$ 1/\sqrt{\epsilon} \f$
 *
 * The bound on the flux coefficient is derived by inverting this value, so it has to be exact.
 */
constexpr double regularInvSqrtPeak(const bool useLinear)
{ return useLinear ? 1.5 : 1.25; }

/*!
 * \ingroup LongWaveModel
 * \brief \f$ 1/\sqrt{x} \f$, continued below the threshold \f$ \epsilon \f$ such that it stays bounded
 *        as \f$ x \to 0 \f$
 *
 * The function is exact above the threshold and C1 across it, so the flux it scales stays
 * monotone and its derivative cannot change sign there. The continuation is linear or
 * quadratic in \f$ x \f$. A non-positive threshold disables the continuation.
 */
template<class Scalar>
Scalar regularInvSqrt(const Scalar x, const Scalar threshold, const bool useLinear)
{
    using std::sqrt;
    if (threshold <= 0.0 || x > threshold)
        return 1.0/sqrt(x);

    const auto f = [](const auto& x){ return 1.0/sqrt(x); };
    const auto dfdx = [](const auto& x){ return -0.5*sqrt(x)/(x*x); };

    if (useLinear)
    {
        const auto a = dfdx(threshold);
        const auto c = f(threshold) - a*threshold;
        return a*x + c;
    }

    const auto a = dfdx(threshold)/(2*threshold);
    const auto c = f(threshold) - a*threshold*threshold;
    return a*x*x + c;
}

/*!
 * \ingroup LongWaveModel
 * \brief Whether `regularInvSqrt` is continued linearly (parameter
 *        `LongWave.Regularization.GradHUseLinear`, default false) or quadratically
 */
inline bool regularizationUsesLinear()
{
    static const bool useLinear = getParam<bool>("LongWave.Regularization.GradHUseLinear", false);
    return useLinear;
}

/*!
 * \ingroup LongWaveModel
 * \brief Lower bound on \f$ |\nabla H| \f$ that bounds the flux coefficient \f$ D \f$ of
 *        \f$ q = -D \nabla H \f$ by `maxDiffusivity`
 *
 * \f$ D = K/(n\sqrt{|\nabla H|}) \f$ with the conveyance \f$ K \f$ diverges as the free surface
 * levels, which is the permanent state of a lake. A fixed lower bound on \f$ |\nabla H| \f$ bounds
 * \f$ D \f$, but since \f$ D \f$ scales with the conveyance, the resulting bound on \f$ D \f$ spans
 * several orders of magnitude between sheet flow and a lake, and it is largest where the water
 * is deepest and most certainly level. Inverting a bound on \f$ D \f$ instead gives a lower bound
 * on the gradient that follows the conveyance: inactive in shallow water, and binding only on
 * deep, near-level water.
 *
 * `peakFactor` is the peak of the regularization of \f$ 1/\sqrt{x} \f$ the bound is used with,
 * relative to \f$ 1/\sqrt{\epsilon} \f$: `regularInvSqrtPeak` for `regularInvSqrt`, and one for
 * the additive shift \f$ 1/\sqrt{x + \epsilon} \f$. An infinite `maxDiffusivity` returns zero.
 */
template<class Scalar>
Scalar diffusivityLimitedThreshold(const Scalar conveyance, const Scalar manningN,
                                   const Scalar maxDiffusivity, const Scalar peakFactor)
{
    if (std::isinf(maxDiffusivity))
        return 0.0;

    const auto root = peakFactor*conveyance/(manningN*maxDiffusivity);
    return root*root;
}

/*!
 * \ingroup LongWaveModel
 * \brief Bound on the flux coefficient in \f$ \mathrm{m^2/s} \f$ (parameter
 *        `LongWave.Regularization.MaxDiffusivity`, default infinite)
 */
template<class Scalar>
Scalar maxDiffusivity()
{
    static const Scalar value = getParam<Scalar>(
        "LongWave.Regularization.MaxDiffusivity", std::numeric_limits<Scalar>::infinity()
    );
    return value;
}

/*!
 * \ingroup LongWaveModel
 * \brief The threshold applied to \f$ |\nabla H| \f$: the more restrictive of the two configured bounds
 *
 * `LongWave.Regularization.GradHEpsilon` (default 1e-8) bounds the free-surface gradient
 * directly, `LongWave.Regularization.MaxDiffusivity` bounds the flux coefficient. Both are
 * lower bounds on the same quantity, so the larger one applies.
 */
template<class Scalar>
Scalar gradHThreshold(const Scalar conveyance, const Scalar manningN)
{
    static const Scalar gradHEps = getParam<Scalar>("LongWave.Regularization.GradHEpsilon", 1e-8);

    using std::max;
    return max(gradHEps, diffusivityLimitedThreshold(
        conveyance, manningN, maxDiffusivity<Scalar>(),
        Scalar(regularInvSqrtPeak(regularizationUsesLinear()))
    ));
}

/*!
 * \ingroup LongWaveModel
 * \brief The depth-dependent part of Manning's law per unit width, \f$ h R^{2/3} \f$ with the
 *        hydraulic radius \f$ R \f$ of unconfined sheet flow
 */
template<class Scalar>
Scalar conveyance(const Scalar h)
{
    using std::max, std::pow;
    return max(Scalar(0.0), h)*pow(regularizedDepth(h), 2.0/3.0);
}

/*!
 * \ingroup LongWaveModel
 * \brief The flux coefficient \f$ D \f$ of \f$ q = -D \nabla H \f$, bounded as the free surface levels
 */
template<class Scalar>
Scalar diffusivity(const Scalar conveyance, const Scalar manningN, const Scalar normGradH)
{
    return conveyance/manningN*regularInvSqrt(
        normGradH, gradHThreshold(conveyance, manningN), regularizationUsesLinear()
    );
}

} // end namespace Dumux::LongWave

#endif
