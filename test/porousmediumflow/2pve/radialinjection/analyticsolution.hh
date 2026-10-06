// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Similarity solution for the injection of a less dense and less viscous fluid into a confined aquifer.
 */

#ifndef DUMUX_TEST_TWOPVE_RADIAL_INJECTION_ANALYTIC_SOLUTION_HH
#define DUMUX_TEST_TWOPVE_RADIAL_INJECTION_ANALYTIC_SOLUTION_HH

#include <cmath>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVETests
 * \brief Similarity solution for radial injection into a confined aquifer with a sharp interface
 *
 * The injected fluid forms a plume of thickness \f$ h(r, t) \f$ below the top of an aquifer of height \f$ H \f$.
 * For a negligible gravity number \f$ \Gamma = 2 \pi \Delta\varrho g k \lambda_w H^2 / Q \f$, the dimensionless
 * plume thickness depends only on the similarity variable \f$ \chi = 2 \pi H \phi (1 - S_{res}) r^2 / (Q t) \f$,
 * \f[
 * \frac{h}{H} =
 * \begin{cases}
 *   1 & \chi \leq 2/\lambda, \\
 *   \frac{1}{\lambda - 1} \left( \sqrt{\frac{2\lambda}{\chi}} - 1 \right) & 2/\lambda < \chi < 2\lambda, \\
 *   0 & \chi \geq 2\lambda,
 * \end{cases}
 * \f]
 * eq. (14) in \cite NordbottenCelia2006, where \f$ \lambda = \lambda_n/\lambda_w \geq 1 \f$ is the mobility ratio,
 * \f$ Q \f$ the volumetric injection rate, \f$ \phi \f$ the porosity, \f$ S_{res} \f$ the residual saturation of the
 * resident fluid, \f$ k \f$ the permeability, \f$ \lambda_w \f$ the mobility of the resident fluid and
 * \f$ \Delta\varrho \f$ the density difference.
 */
template<class Scalar>
class TwoPVERadialInjectionSimilaritySolution
{
public:
    /*!
     * \param mobilityRatio        mobility of the injected fluid divided by the mobility of the resident fluid
     * \param aquiferHeight        height of the aquifer
     * \param porosity             porosity of the aquifer
     * \param residualSaturation   residual saturation of the resident fluid in the plume
     * \param injectionRate        volumetric injection rate
     */
    TwoPVERadialInjectionSimilaritySolution(Scalar mobilityRatio,
                                            Scalar aquiferHeight,
                                            Scalar porosity,
                                            Scalar residualSaturation,
                                            Scalar injectionRate)
    : mobilityRatio_(mobilityRatio)
    , aquiferHeight_(aquiferHeight)
    , porosity_(porosity)
    , residualSaturation_(residualSaturation)
    , injectionRate_(injectionRate)
    {
        if (mobilityRatio_ < 1.0)
            DUNE_THROW(Dune::InvalidStateException, "The similarity solution requires a mobility ratio of at least one");
    }

    //! The similarity variable \f$ \chi \f$ at radius r and time t
    Scalar similarityVariable(Scalar radius, Scalar time) const
    {
        return 2.0*M_PI*aquiferHeight_*porosity_*(1.0 - residualSaturation_)*radius*radius/(injectionRate_*time);
    }

    //! The dimensionless plume thickness \f$ h/H \f$ as a function of the similarity variable
    Scalar dimensionlessPlumeThickness(Scalar chi) const
    {
        using std::sqrt;
        if (chi <= 2.0/mobilityRatio_)
            return 1.0;
        else if (chi >= 2.0*mobilityRatio_)
            return 0.0;
        else
            return (sqrt(2.0*mobilityRatio_/chi) - 1.0)/(mobilityRatio_ - 1.0);
    }

    //! The gas plume distance above the bottom of the aquifer at radius r and time t
    Scalar gasPlumeDistance(Scalar radius, Scalar time) const
    {
        return aquiferHeight_*(1.0 - dimensionlessPlumeThickness(similarityVariable(radius, time)));
    }

    //! The radius of the plume tip at time t
    Scalar plumeExtent(Scalar time) const
    {
        using std::sqrt;
        return sqrt(2.0*mobilityRatio_*injectionRate_*time/(2.0*M_PI*aquiferHeight_*porosity_*(1.0 - residualSaturation_)));
    }

private:
    Scalar mobilityRatio_;
    Scalar aquiferHeight_;
    Scalar porosity_;
    Scalar residualSaturation_;
    Scalar injectionRate_;
};

} // end namespace Dumux

#endif
