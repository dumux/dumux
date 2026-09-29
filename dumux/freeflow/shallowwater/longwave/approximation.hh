// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \brief The long-wave approximations of the shallow water equations and their range of validity
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_APPROXIMATION_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_APPROXIMATION_HH

#include <algorithm>
#include <cmath>
#include <string>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux::LongWave {

/*!
 * \ingroup LongWaveModel
 * \brief Which long-wave approximation of the shallow water equations to solve
 *
 * Selected by the parameter `LongWave.Approximation`. The flux is driven by
 * \f$ \nabla z + w \nabla h \f$, so the approximations differ only in the weight \f$ w \f$.
 */
enum class WaveApproximation
{
    diffusive, //!< \f$ w = 1 \f$: driven by the free-surface gradient
    kinematic, //!< \f$ w = 0 \f$: driven by the bed slope alone
    inertiaCorrected //!< \f$ w = \max(0, 1 - \mathsf{V}^2) \f$ with the Vedernikov number \f$ \mathsf{V} \f$
};

/*!
 * \ingroup LongWaveModel
 * \brief The approximation selected by the parameter `LongWave.Approximation`
 *        (`diffusive` (default), `kinematic` or `inertiacorrected`)
 */
inline WaveApproximation waveApproximation()
{
    static const auto approximation = []
    {
        const auto name = getParam<std::string>("LongWave.Approximation", "diffusive");
        if (name == "diffusive") return WaveApproximation::diffusive;
        if (name == "kinematic") return WaveApproximation::kinematic;
        if (name == "inertiacorrected") return WaveApproximation::inertiaCorrected;
        DUNE_THROW(ParameterException, "Unknown LongWave.Approximation '" << name
                   << "', expected 'diffusive', 'kinematic' or 'inertiacorrected'");
    }();
    return approximation;
}

/*!
 * \ingroup LongWaveModel
 * \brief Froude number of normal flow under Manning friction in a wide channel, with
 *        \f$ g = 9.81\,\mathrm{m/s^2} \f$
 *
 * Taking the Froude number for normal flow on the local bed slope, rather than on the free
 * surface, keeps it independent of the free-surface gradient, so it stays defined as the
 * surface levels. Since \f$ \mathsf{Fr} \propto h^{1/6} \f$, the dependence on the depth is weak.
 */
template<class Scalar>
Scalar froudeNumber(const Scalar bedSlope, const Scalar manningN, const Scalar depth)
{
    using std::sqrt, std::pow, std::max;
    static constexpr Scalar gravity = 9.81;
    return sqrt(bedSlope)/(manningN*sqrt(gravity))*pow(max(depth, Scalar(0.0)), 1.0/6.0);
}

/*!
 * \ingroup LongWaveModel
 * \brief Vedernikov number of normal flow, \f$ \mathsf{V} = \frac{2}{3}\mathsf{Fr} \f$ for Manning
 *        friction in a wide channel
 *
 * Above one, the kinematic wave travels faster than the fastest gravity wave and the flow is
 * past the roll-wave threshold. There, the smoothing applied by the diffusive wave is spurious,
 * which makes this number the criterion for whether the long-wave approximations apply at all.
 */
template<class Scalar>
Scalar vedernikovNumber(const Scalar bedSlope, const Scalar manningN, const Scalar depth)
{ return 2.0/3.0*froudeNumber(bedSlope, manningN, depth); }

/*!
 * \ingroup LongWaveModel
 * \brief Vedernikov number of normal flow computed from a discharge per unit width
 *
 * The discharge enters as \f$ q^{1/10} \f$, so the number is set almost entirely by the bed
 * slope and the roughness. A rough estimate of the discharge suffices to map it over a
 * catchment before anything is solved.
 */
template<class Scalar>
Scalar vedernikovNumberFromDischarge(const Scalar bedSlope, const Scalar manningN,
                                     const Scalar dischargePerWidth)
{
    using std::sqrt, std::pow, std::max;
    const auto depth = pow(max(dischargePerWidth, Scalar(0.0))*manningN/sqrt(bedSlope), 3.0/5.0);
    return vedernikovNumber(bedSlope, manningN, depth);
}

/*!
 * \ingroup LongWaveModel
 * \brief Weight \f$ w \f$ of the depth gradient in the flux driver \f$ \nabla z + w \nabla h \f$
 *
 * Linearizing the shallow water equations about uniform flow gives the diffusivity
 * \f$ q/(2 S_0) (1 - \mathsf{V}^2) \f$. The diffusive wave keeps the free-surface gradient but drops
 * the inertial terms that supply the factor \f$ 1 - \mathsf{V}^2 \f$, so it over-diffuses as
 * \f$ \mathsf{V} \f$ approaches one. The inertia-corrected approximation restores the factor,
 * which turns the diffusive wave into the kinematic wave where the latter is the better
 * approximation.
 *
 * The weight is clamped at zero since a negative weight makes the equation backward-parabolic.
 * Beyond \f$ \mathsf{V} = 1 \f$ the physical flow develops roll waves, which no long-wave
 * approximation reproduces.
 */
template<class Scalar>
Scalar freeSurfaceWeight(const Scalar bedSlope, const Scalar depth, const Scalar manningN)
{
    switch (waveApproximation())
    {
        case WaveApproximation::diffusive: return 1.0;
        case WaveApproximation::kinematic: return 0.0;
        case WaveApproximation::inertiaCorrected: break;
    }

    using std::max;
    const auto vedernikov = vedernikovNumber(bedSlope, manningN, depth);
    return max(Scalar(0.0), 1.0 - vedernikov*vedernikov);
}

} // end namespace Dumux::LongWave

#endif
