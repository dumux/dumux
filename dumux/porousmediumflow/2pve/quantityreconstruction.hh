// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief Reconstructs quantities from the coarse to the fine level of the VE scheme.
 */

#ifndef DUMUX_TWOPVE_QUANTITY_RECONSTRUCTION_HH
#define DUMUX_TWOPVE_QUANTITY_RECONSTRUCTION_HH

#include <algorithm>
#include <array>
#include <cmath>

#include <dune/common/exceptions.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/nonlinear/findscalarroot.hh>

namespace Dumux {

namespace TwoPVE {

/*!
 * \ingroup TwoPVEModel
 * \brief The gas plume distance of a column and its minimum over all previous time steps
 */
template<typename Scalar>
struct GasPlumeDistancesData
{
    Scalar current;
    Scalar minimum;
};

/*!
 * \ingroup TwoPVEModel
 * \brief The densities of the wetting and the nonwetting phase
 */
template<typename Scalar>
struct PhaseDensitiesData
{
    Scalar wetting;
    Scalar nonwetting;
};

/*!
 * \ingroup TwoPVEModel
 * \brief The viscosities of the wetting and the nonwetting phase
 */
template<typename Scalar>
struct PhaseViscositiesData
{
    Scalar wetting;
    Scalar nonwetting;
};

/*!
 * \ingroup TwoPVEModel
 * \brief The residual saturations of the wetting and the nonwetting phase
 */
template<typename Scalar>
struct ResidualSaturationsData
{
    Scalar wetting;
    Scalar nonwetting;
};

/*!
 * \ingroup TwoPVEModel
 * \brief The parameters of the Brooks-Corey material law
 */
template<typename Scalar>
struct BrooksCoreyParametersData
{
    Scalar lambda;
    Scalar entryPressure;
};

} // end namespace TwoPVE

/*!
 * \ingroup TwoPVEModel
 * \brief Reconstructs the fine-level quantities of a column from its coarse-level quantities assuming vertical equilibrium
 *
 * The reconstruction uses the Brooks-Corey material law \cite brooks1964hydrau.
 * Below the gas plume distance, the pore space is water-saturated apart from residually trapped gas.
 * Above the gas plume distance, the saturation follows from the capillary pressure in hydrostatic equilibrium.
 * The averages of the saturation and the relative permeabilities over fine-level cells are evaluated in closed form.
 *
 * \tparam Scalar the scalar type
 * \tparam FluidSystem the immiscible two-phase fluid system
 */
template<class Scalar, class FluidSystem>
class TwoPVEQuantityReconstruction
{
    static constexpr int wettingPhaseIdx = FluidSystem::phase0Idx;
    static constexpr int nonwettingPhaseIdx = FluidSystem::phase1Idx;
    static constexpr int numPhases = FluidSystem::numPhases;

    using GasPlumeDistances = TwoPVE::GasPlumeDistancesData<Scalar>;
    using PhaseDensities = TwoPVE::PhaseDensitiesData<Scalar>;
    using PhaseViscosities = TwoPVE::PhaseViscositiesData<Scalar>;
    using ResidualSaturations = TwoPVE::ResidualSaturationsData<Scalar>;
    using BrooksCoreyParameters = TwoPVE::BrooksCoreyParametersData<Scalar>;

public:

    /*!
     * \brief Computes the gas plume distance for a coarse column, which is the height of the gas plume relative to the bottom of the domain. The gas plume distance should be a value in [0, domainHeight]
     *
     * The gas plume distance is the root of the column water balance. The water content of the column
     * accounts for the trapped-gas region between the minimum gas plume distance and the gas plume distance,
     * such that the column integral of the reconstructed saturation equals the coarse-level saturation.
     * If the coarse-level saturation exceeds the water content of a column without mobile gas,
     * the gas plume distance is the domain height.
     *
     * \param densities             contains the phase densities (here: 2 phases)
     * \param residualSaturations   contains the phase residual saturation (here: 2 phases)
     * \param gravityNorm           norm of the gravity
     * \param domainHeight          height of the whole domain
     * \param saturationWCoarse     wetting-phase saturation (on the coarse level)
     * \param minGasPlumeDist       minimum gas plume distance of the previous time steps, a value in [0, domainHeight]
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     * \throws NumericalProblem if no gas plume distance in [0, domainHeight] matches the coarse-level saturation
     */
    Scalar computeGasPlumeDist(const PhaseDensities& densities,
                               const ResidualSaturations& residualSaturations,
                               Scalar gravityNorm,
                               Scalar domainHeight,
                               Scalar saturationWCoarse,
                               Scalar minGasPlumeDist,
                               const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar densityW = densities.wetting;
        const Scalar densityNw = densities.nonwetting;
        const Scalar snr = residualSaturations.nonwetting;

        if (gravityNorm == 0.0)
            DUNE_THROW(Dune::InvalidStateException, "The two-phase vertical-equilibrium model requires nonzero gravity because its reconstruction assumes gravity-driven vertical segregation");

        if (densityW <= densityNw)
            DUNE_THROW(Dune::InvalidStateException, "The two-phase VE reconstruction requires the wetting phase to be denser than the nonwetting phase. Got rho_w=" << densityW << " and rho_n=" << densityNw);

        // water content of the column minus the coarse-level water content, monotonically increasing in the gas plume distance
        using std::max; using std::min;
        const auto massConservation = [&](const Scalar gasPlumeDist)
        {
            const Scalar waterBelowMinimum = min(gasPlumeDist, minGasPlumeDist);
            const Scalar waterTrappedRegion = (1.0-snr) * max(gasPlumeDist - minGasPlumeDist, 0.0);
            const Scalar waterAbovePlume = integrateSaturationWAbovePlume_(gasPlumeDist, domainHeight, gasPlumeDist, densities, residualSaturations, gravityNorm, brooksCoreyParameters);
            return waterBelowMinimum + waterTrappedRegion + waterAbovePlume - saturationWCoarse * domainHeight;
        };

        if (massConservation(domainHeight) <= 0.0)
            return domainHeight;

        const Scalar residualAtBottom = massConservation(0.0);
        if (residualAtBottom == 0.0)
            return 0.0;
        if (residualAtBottom > 0.0)
            DUNE_THROW(NumericalProblem, "No gas plume distance in [0, " << domainHeight << "] matches the coarse-level wetting-phase saturation " << saturationWCoarse
                                         << " (densities: " << densityW << ", " << densityNw << ")");

        return findScalarRootBrent(0.0, domainHeight, massConservation);
    }

    /*!
     * \brief Computes the capillary pressure on the coarse level
     *
     * \param gasPlumeDist    value of the gas plume distance
     * \param densities       contains the phase densities (here: 2 phases)
     * \param gravityNorm     norm of the gravity
     * \param entryPressureBC entry pressure of the Brooks-Corey model
     */
    Scalar computeCapillaryPressureCoarse(Scalar gasPlumeDist,
                                          const PhaseDensities& densities,
                                          Scalar gravityNorm,
                                          Scalar entryPressureBC) const
    {
        return gravityNorm * gasPlumeDist * (densities.nonwetting - densities.wetting) + entryPressureBC;
    }

    /*!
     * \brief Reconstructs the fine-level wetting-phase and nonwetting-phase pressures
     *
     * \param gasPlumeDist      value of the gas plume distance
     * \param densities         contains the phase densities (here: 2 phases)
     * \param gravityNorm       norm of the gravity
     * \param heightAboveBottom height (relative to the bottom of the domain) at which the pressures should be reconstructed
     * \param pressureWCoarse   wetting-phase pressure on the coarse level
     * \param entryPressureBC   entry pressure of the Brooks-Corey model
     */
    std::array<Scalar, numPhases> reconstructPressure(Scalar gasPlumeDist,
                                                      const PhaseDensities& densities,
                                                      Scalar gravityNorm,
                                                      Scalar heightAboveBottom,
                                                      Scalar pressureWCoarse,
                                                      Scalar entryPressureBC) const
    {
        const Scalar referenceDensityW = densities.wetting;
        const Scalar referenceDensityNw = densities.nonwetting;

        // the nonwetting-phase pressure below the gas plume exceeds the wetting-phase pressure by the entry pressure (capillary fringe model)
        std::array<Scalar, numPhases> reconstructedPressures;
        reconstructedPressures[wettingPhaseIdx] = pressureWCoarse - referenceDensityW * gravityNorm * heightAboveBottom;
        if (heightAboveBottom <= gasPlumeDist)
            reconstructedPressures[nonwettingPhaseIdx] = reconstructedPressures[wettingPhaseIdx] + entryPressureBC;
        else
            reconstructedPressures[nonwettingPhaseIdx] = pressureWCoarse - referenceDensityW * gravityNorm * gasPlumeDist - referenceDensityNw * gravityNorm * (heightAboveBottom - gasPlumeDist) + entryPressureBC;

        return reconstructedPressures;
    }

    /*!
     * \brief Reconstructs the fine-level wetting-phase saturation as a mean integral of the saturation across the cell height
     *
     * Compute the saturation as the mean of the reconstructed saturation over the height of a fine-level cell.
     *
     * \param gasPlumeDistances     contains the gas plume distance and the minimum gas plume distance
     * \param densities             contains the phase densities (here: 2 phases)
     * \param residualSaturations   contains the phase residual saturation (here: 2 phases)
     * \param gravityNorm           norm of the gravity
     * \param heightAboveBottom     height (relative to the bottom of the domain) at which the saturation should be reconstructed
     * \param deltaZ                discretization width of the fine-level grid in vertical direction
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     */
    Scalar reconstructSaturation(const GasPlumeDistances& gasPlumeDistances,
                                 const PhaseDensities& densities,
                                 const ResidualSaturations& residualSaturations,
                                 Scalar gravityNorm,
                                 Scalar heightAboveBottom,
                                 Scalar deltaZ,
                                 const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar gasPlumeDist = gasPlumeDistances.current;
        const Scalar minGasPlumeDist = gasPlumeDistances.minimum;

        const Scalar lowerBound = heightAboveBottom - deltaZ/2.0;
        const Scalar upperBound = heightAboveBottom + deltaZ/2.0;
        const Scalar cellHeight = upperBound - lowerBound;

        const Scalar snr = residualSaturations.nonwetting;
        const Scalar saturationWBelowPlume = 1.0;
        const Scalar saturationWTrappedRegion = 1.0 - snr;
        const auto integrateSaturationWAbovePlume = [&](const Scalar lower, const Scalar upper)
        {
            return integrateSaturationWAbovePlume_(lower, upper, gasPlumeDist, densities, residualSaturations, gravityNorm, brooksCoreyParameters);
        };

        Scalar saturationIntegral = 0.0;

        // differentiate between six cases:
        //     1) fine-level element is completely under the gas plume,
        //     2) fine-level element is completely above the gas plume,
        //     3) fine-level element is completely in trapped gas region
        //     4) fine-level element intersects only historical minimum plume distance
        //     5) fine-level element intersects only the current gas plume distance
        //     6) fine-level element intersects both plume distances

        if (upperBound <= minGasPlumeDist) // 1)
            saturationIntegral = saturationWBelowPlume*cellHeight;
        else if (lowerBound >= gasPlumeDist) // 2)
            saturationIntegral = integrateSaturationWAbovePlume(lowerBound, upperBound);
        else if (lowerBound >= minGasPlumeDist && upperBound <= gasPlumeDist) // 3)
            saturationIntegral = saturationWTrappedRegion*cellHeight;
        else if (lowerBound < minGasPlumeDist && upperBound <= gasPlumeDist) // 4)
        {
            const Scalar integralBelowMinimum = saturationWBelowPlume * (minGasPlumeDist - lowerBound);
            const Scalar integralTrappedRegion = saturationWTrappedRegion * (upperBound - minGasPlumeDist);
            saturationIntegral = integralBelowMinimum + integralTrappedRegion;
        }
        else if (lowerBound >= minGasPlumeDist && upperBound > gasPlumeDist) // 5)
        {
            const Scalar integralTrappedRegion = saturationWTrappedRegion * (gasPlumeDist - lowerBound);
            const Scalar integralAbovePlume = integrateSaturationWAbovePlume(gasPlumeDist, upperBound);
            saturationIntegral = integralTrappedRegion + integralAbovePlume;
        }
        else // 6)
        {
            const Scalar integralBelowMinimum = saturationWBelowPlume * (minGasPlumeDist - lowerBound);
            const Scalar integralTrappedRegion = saturationWTrappedRegion * (gasPlumeDist - minGasPlumeDist);
            const Scalar integralAbovePlume = integrateSaturationWAbovePlume(gasPlumeDist, upperBound);
            saturationIntegral = integralBelowMinimum + integralTrappedRegion + integralAbovePlume;
        }

        return saturationIntegral/cellHeight;
    }

    /*!
     * \brief Reconstructs the fine-level capillary pressure
     *
     * \param gasPlumeDist      value of the gas plume distance
     * \param densities         contains the phase densities (here: 2 phases)
     * \param gravityNorm       norm of the gravity
     * \param heightAboveBottom height (relative to the bottom of the domain) at which the capillary pressure should be reconstructed
     * \param entryPressureBC   entry pressure of the Brooks-Corey model
     */
    Scalar reconstructCapillaryPressure(Scalar gasPlumeDist,
                                        const PhaseDensities& densities,
                                        Scalar gravityNorm,
                                        Scalar heightAboveBottom,
                                        Scalar entryPressureBC) const
    {
        const Scalar referenceDensityW = densities.wetting;
        const Scalar referenceDensityNw = densities.nonwetting;

        if (heightAboveBottom <= gasPlumeDist)
            return entryPressureBC;
        else
            return referenceDensityW * gravityNorm * (heightAboveBottom - gasPlumeDist) + entryPressureBC - referenceDensityNw * gravityNorm * (heightAboveBottom - gasPlumeDist);
    }

    /*!
     * \brief Reconstructs the fine-level wetting-phase and nonwetting-phase mobilities
     *
     * \param gasPlumeDistances     contains the gas plume distance and the minimum gas plume distance
     * \param densities             contains the phase densities (here: 2 phases)
     * \param viscosities           contains the phase viscosities (here: 2 phases)
     * \param gravityNorm           norm of the gravity
     * \param heightAboveBottom     height (relative to the bottom of the domain) at which the mobilities should be reconstructed
     * \param deltaZ                discretization width of the fine-level grid in vertical direction
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     */
    std::array<Scalar, numPhases> reconstructMobilitiesFine(const GasPlumeDistances& gasPlumeDistances,
                                                            const PhaseDensities& densities,
                                                            const PhaseViscosities& viscosities,
                                                            Scalar gravityNorm,
                                                            Scalar heightAboveBottom,
                                                            Scalar deltaZ,
                                                            const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar gasPlumeDist = gasPlumeDistances.current;
        const Scalar lowerBound = heightAboveBottom - deltaZ/2.0;
        const Scalar upperBound = heightAboveBottom + deltaZ/2.0;

        const Scalar relPermWBelowPlume = 1.0;
        const Scalar relPermNwBelowPlume = 0.0;

        // Brooks-Corey relative permeabilities in terms of the effective saturation Se = u^(-lambda):
        // krw = Se^(2/lambda + 3) and krn = (1 - Se)^2 (1 - Se^(2/lambda + 1)), expanded into powers of u
        const Scalar lambdaBC = brooksCoreyParameters.lambda;
        const auto integratePower = [&](const Scalar exponent, const Scalar lower, const Scalar upper)
        {
            return integratePowerAbovePlume_(exponent, lower, upper, gasPlumeDist, densities, gravityNorm, brooksCoreyParameters.entryPressure);
        };
        const auto integrateRelPermWAbovePlume = [&](const Scalar lower, const Scalar upper)
        {
            return integratePower(2.0 + 3.0*lambdaBC, lower, upper);
        };
        const auto integrateRelPermNwAbovePlume = [&](const Scalar lower, const Scalar upper)
        {
            const Scalar integral = (upper - lower)
                                    - 2.0*integratePower(lambdaBC, lower, upper)
                                    + integratePower(2.0*lambdaBC, lower, upper)
                                    - integratePower(2.0 + lambdaBC, lower, upper)
                                    + 2.0*integratePower(2.0 + 2.0*lambdaBC, lower, upper)
                                    - integratePower(2.0 + 3.0*lambdaBC, lower, upper);

            // close to the gas plume distance, where krn vanishes, the terms cancel up to round-off errors of either sign
            using std::max;
            return max(integral, 0.0);
        };

        // the relative permeabilities are piecewise defined, with a kink or jump at the gas plume distance
        const auto integrateRelativePermeability = [&](const auto& integralAbovePlume, const Scalar relPermBelowPlume)
        {
            const Scalar cellHeight = upperBound - lowerBound;

            // differentiate between three cases:
            //     1) fine-level element is completely under the gas plume,
            //     2) fine-level element is completely above the gas plume,
            //     3) fine-level element intersects with gas plume
            //     for zp_min<z<zp, there might be snr present and the rel perm for gas is treated as 0, thus the check z<zp is enough
            if (upperBound <= gasPlumeDist)
                return relPermBelowPlume * cellHeight;
            else if (lowerBound >= gasPlumeDist)
                return integralAbovePlume(lowerBound, upperBound);
            else
                return relPermBelowPlume * (gasPlumeDist - lowerBound) + integralAbovePlume(gasPlumeDist, upperBound);
        };

        // cell averages of the relative permeabilities divided by the viscosities
        std::array<Scalar, numPhases> mobilities;
        mobilities[wettingPhaseIdx] = integrateRelativePermeability(integrateRelPermWAbovePlume, relPermWBelowPlume)/(upperBound - lowerBound)/viscosities.wetting;
        mobilities[nonwettingPhaseIdx] = integrateRelativePermeability(integrateRelPermNwAbovePlume, relPermNwBelowPlume)/(upperBound - lowerBound)/viscosities.nonwetting;
        return mobilities;
    }

private:
    /*!
     * \brief Integrates \f$ u^{-m} \f$ with \f$ u = 1 + \Delta\varrho g (z - z_p)/p_e \f$ over \f$ z \in [z_l, z_u] \f$ above the gas plume distance \f$ z_p \f$
     *
     * Above the gas plume distance, the Brooks-Corey effective wetting-phase saturation is \f$ u^{-\lambda} \f$,
     * where \f$ \Delta\varrho \f$ is the density difference of the phases, \f$ g \f$ the norm of the gravity
     * and \f$ p_e \f$ the entry pressure.
     *
     * \param exponent        the exponent \f$ m \f$
     * \param lowerBound      lower integration bound \f$ z_l \geq z_p \f$
     * \param upperBound      upper integration bound \f$ z_u \geq z_l \f$
     * \param gasPlumeDist    gas plume distance \f$ z_p \f$
     * \param densities       contains the phase densities (here: 2 phases)
     * \param gravityNorm     norm of the gravity
     * \param entryPressureBC entry pressure of the Brooks-Corey model
     */
    Scalar integratePowerAbovePlume_(Scalar exponent,
                                     Scalar lowerBound,
                                     Scalar upperBound,
                                     Scalar gasPlumeDist,
                                     const PhaseDensities& densities,
                                     Scalar gravityNorm,
                                     Scalar entryPressureBC) const
    {
        using std::expm1; using std::log1p; using std::pow;
        const Scalar lengthScale = entryPressureBC/((densities.wetting - densities.nonwetting)*gravityNorm);
        const Scalar uLower = 1.0 + (lowerBound - gasPlumeDist)/lengthScale;
        const Scalar logRatio = log1p((upperBound - lowerBound)/(lengthScale*uLower));
        const Scalar k = 1.0 - exponent;
        if (k == 0.0)
            return lengthScale*logRatio;

        // expm1 keeps the antiderivative accurate for exponents close to one, where it approaches the logarithm
        return lengthScale*pow(uLower, k)*expm1(k*logRatio)/k;
    }

    /*!
     * \brief Integrates the wetting-phase saturation over \f$ [z_l, z_u] \f$ above the gas plume distance
     *
     * \param lowerBound            lower integration bound, not below the gas plume distance
     * \param upperBound            upper integration bound
     * \param gasPlumeDist          gas plume distance
     * \param densities             contains the phase densities (here: 2 phases)
     * \param residualSaturations   contains the phase residual saturation (here: 2 phases)
     * \param gravityNorm           norm of the gravity
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     */
    Scalar integrateSaturationWAbovePlume_(Scalar lowerBound,
                                           Scalar upperBound,
                                           Scalar gasPlumeDist,
                                           const PhaseDensities& densities,
                                           const ResidualSaturations& residualSaturations,
                                           Scalar gravityNorm,
                                           const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar swr = residualSaturations.wetting;
        const Scalar snr = residualSaturations.nonwetting;
        const Scalar integralEffectiveSaturation = integratePowerAbovePlume_(brooksCoreyParameters.lambda, lowerBound, upperBound, gasPlumeDist,
                                                                             densities, gravityNorm, brooksCoreyParameters.entryPressure);
        return swr*(upperBound - lowerBound) + (1.0 - swr - snr)*integralEffectiveSaturation;
    }
};

} // end namespace Dumux

#endif
