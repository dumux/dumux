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
 * \brief Froude number of normal flow under Manning friction in a wide channel,
 *        \f$ \mathsf{Fr} = \sqrt{S_0}\, h^{1/6}/(n \sqrt{g}) \f$ with \f$ g = 9.81\,\mathrm{m/s^2} \f$
 *
 * Evaluated with the bed slope \f$ S_0 \f$, so that it does not depend on the free-surface gradient.
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
 * Above one, the kinematic wave is faster than the downstream gravity wave and uniform flow is
 * unstable to roll waves, which the long-wave approximations do not represent.
 */
template<class Scalar>
Scalar vedernikovNumber(const Scalar bedSlope, const Scalar manningN, const Scalar depth)
{ return 2.0/3.0*froudeNumber(bedSlope, manningN, depth); }

/*!
 * \ingroup LongWaveModel
 * \brief Vedernikov number of normal flow computed from a discharge per unit width
 *
 * The normal depth follows from Manning's law, \f$ h = (q n/\sqrt{S_0})^{3/5} \f$.
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
 * Linearized about uniform flow with the discharge per unit width \f$ q \f$ and the bed slope
 * \f$ S_0 \f$, the shallow water equations have the diffusivity \f$ q/(2 S_0) (1 - \mathsf{V}^2) \f$
 * and the diffusive wave has \f$ q/(2 S_0) \f$. The inertia-corrected weight recovers the former.
 * It is clamped at zero, since a negative weight makes the equation backward-parabolic.
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
