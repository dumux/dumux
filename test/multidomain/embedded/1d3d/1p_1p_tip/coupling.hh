// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Coupling coefficients of the vessel tip benchmark for the different coupling modes.
 */
#ifndef DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_COUPLING_HH
#define DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_COUPLING_HH

#include <cmath>
#include <map>
#include <vector>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

#include "exactsolution.hh"

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The source distribution in the tissue that corresponds to a 1d-3d coupling mode
 */
template<class CouplingManager>
constexpr VesselSourceDistribution vesselSourceDistribution()
{
    if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::kernel)
        return VesselSourceDistribution::kernel;
    else if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::surface)
        return VesselSourceDistribution::surface;
    else if constexpr (CouplingManager::couplingMode == Embedded1d3dCouplingMode::average)
        return VesselSourceDistribution::line;
    else
        static_assert(AlwaysFalse<CouplingManager>::value,
                      "The vessel tip benchmark has an exact solution for the average, surface and kernel coupling modes");
}

/*!
 * \ingroup EmbeddedTests
 * \brief The flux scaling factor \f$ \xi \f$ of the kernel method with a constant kernel of radius
 *        \f$ \varrho \f$ for a distance r between a vessel segment and a tissue cell
 *
 * Koch et al. (2020), Eqs. (3.18) and (3.20), for the coupling coefficient \f$ \beta \f$ (unit tissue
 * permeability and viscosity). For r = 0 it is the continuous flux scaling factor.
 */
template<class Scalar>
Scalar kernelFluxScalingFactor(Scalar r, Scalar radius, Scalar kernelRadius, Scalar couplingCoefficient)
{
    using std::log;
    const Scalar shape = r < kernelRadius ? r*r/(2.0*kernelRadius*kernelRadius) + log(kernelRadius/radius) - 0.5
                                          : log(r/radius);
    return 1.0/(1.0 + couplingCoefficient/(2.0*M_PI)*shape);
}

/*!
 * \ingroup EmbeddedTests
 * \brief The coupling coefficient \f$ \beta \f$ for each coupling source id
 *
 * The coefficient is evaluated at the height of the source on the vessel center-line.
 */
template<class CouplingManager, class Scalar>
std::vector<Scalar> vesselTipCouplingCoefficients(const CouplingManager& couplingManager,
                                                  const VesselTipExactSolution<Scalar>& exact)
{
    constexpr auto distribution = vesselSourceDistribution<CouplingManager>();

    std::map<Scalar, Scalar> coefficientAtHeight;
    std::vector<Scalar> coefficients;
    for (const auto& source : couplingManager.lowDimPointSources())
    {
        const Scalar z = source.position()[2];
        auto it = coefficientAtHeight.find(z);
        if (it == coefficientAtHeight.end())
            it = coefficientAtHeight.emplace(z, exact.couplingCoefficient(distribution, z)).first;

        if (source.id() >= coefficients.size())
            coefficients.resize(source.id() + 1, 0.0);
        coefficients[source.id()] = it->second;
    }
    return coefficients;
}

} // end namespace Dumux

#endif
