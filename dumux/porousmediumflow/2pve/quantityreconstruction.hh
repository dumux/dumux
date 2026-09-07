// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVE
 * \brief Reconstructs quantities from the coarse to the fine level of the VE scheme.
 */

#ifndef DUMUX_TWOPVE_QUANTITY_RECONSTRUCTION_HH
#define DUMUX_TWOPVE_QUANTITY_RECONSTRUCTION_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/nonlinear/findscalarroot.hh>
#include <dumux/common/integrate.hh>

namespace Dumux {

namespace TwoPVE {
    template<typename Scalar>
    struct GasPlumeDistancesData
    {
        Scalar current;
        Scalar minimum;
    };

    template<typename Scalar>
    struct PhaseDensitiesData
    {
        Scalar wetting;
        Scalar nonwetting;
    };

    template<typename Scalar>
    struct PhaseViscositiesData
    {
        Scalar wetting;
        Scalar nonwetting;
    };

    template<typename Scalar>
    struct ResidualSaturationsData
    {
        Scalar wetting;
        Scalar nonwetting;
    };

    template<typename Scalar>
    struct BrooksCoreyParametersData
    {
        Scalar lambda;
        Scalar entryPressure;
    };
} // end namespace TwoPVE

template<class TypeTag>
class TwoPVEQuantityReconst
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum {
        waterPhaseIdx = FluidSystem::phase0Idx, // = 0
        gasPhaseIdx = FluidSystem::phase1Idx,   // = 1
        numPhases = FluidSystem::numPhases
    };
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GasPlumeDistances = TwoPVE::GasPlumeDistancesData<Scalar>;
    using PhaseDensities = TwoPVE::PhaseDensitiesData<Scalar>;
    using PhaseViscosities = TwoPVE::PhaseViscositiesData<Scalar>;
    using ResidualSaturations = TwoPVE::ResidualSaturationsData<Scalar>;
    using BrooksCoreyParameters = TwoPVE::BrooksCoreyParametersData<Scalar>;

public:

    TwoPVEQuantityReconst(const std::string& spatialParamsGroup = "")
    {}


    /*!
     * \brief Computes the gas plume distance for a coarse column, which is the height of the gas plume relative to the bottom of the domain. The gas plume distance should be a value in [0, domainHeight]
     *
     * \param densities             contains the phase densities (here: 2 phases)
     * \param residualSaturations   contains the phase residual saturation (here: 2 phases)
     * \param gravityNorm           norm of the gravity
     * \param domainHeight          height of the whole domain
     * \param satWCoarse            wetting-phase saturation (on the coarse level)
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     */
    Scalar computeGasPlumeDist(const PhaseDensities& densities,
                               const ResidualSaturations& residualSaturations,
                               const Scalar& gravityNorm,
                               const Scalar& domainHeight,
                               const Scalar& satWCoarse,
                               const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar densityW = densities.wetting;
        const Scalar densityNw = densities.nonwetting;
        const Scalar swr = residualSaturations.wetting;
        const Scalar snr = residualSaturations.nonwetting;
        const Scalar lambdaBC = brooksCoreyParameters.lambda;
        const Scalar entryPressureBC = brooksCoreyParameters.entryPressure;

        if (float_equal_(gravityNorm, 0.0))
            DUNE_THROW(Dune::InvalidStateException, "The two-phase vertical-equilibrium model requires nonzero gravity because its reconstruction assumes gravity-driven vertical segregation");

        if (float_equal_(densityW - densityNw, 0.0) || densityW<densityNw)
            DUNE_THROW(Dune::InvalidStateException, "The two-phase VE reconstruction requires the wetting phase to be denser than the nonwetting phase. Got rho_w=" << densityW << " and rho_n=" << densityNw);

        if (float_equal_(lambdaBC, 1.0))
            DUNE_THROW(Dune::InvalidStateException, "Brooks-Corey lambda=1 is not supported by the analytical gas-plume-distance formula");

        // lambda function for mass content in column, to be solved
        const auto massConservation = [&](const Scalar gasPlumeDist)
        {
            auto A = std::pow(entryPressureBC,lambdaBC) * (1.0-swr-snr);
            return
            gasPlumeDist - 0.0
            - satWCoarse * domainHeight
            + (1.0/(1.0-lambdaBC)) * (1.0/((densityW - densityNw)*gravityNorm)) * A * ( std::pow(entryPressureBC+(densityW - densityNw)*gravityNorm*(domainHeight-gasPlumeDist), 1.0-lambdaBC) - std::pow(entryPressureBC, 1.0-lambdaBC) )
            + swr*(domainHeight-gasPlumeDist);
        };

        // lambda function for derivative of mass content in column
        const auto massConservationDerivative = [&](const Scalar gasPlumeDist)
        {
            auto A = std::pow(entryPressureBC,lambdaBC) * (1.0-swr-snr);
            return
            1.0
            + A * ( std::pow(entryPressureBC+(densityW - densityNw)*gravityNorm*(domainHeight-gasPlumeDist), -lambdaBC) ) *(-1.0)
            -swr;
        };

        Scalar gasPlumeDistance = domainHeight;
        const Scalar initialGuess = 0.5 * domainHeight; // using initialGuess=domainHeight leads to convergence issues

        if (satWCoarse < 1.0)
        {
            try
            {
                gasPlumeDistance = findScalarRootNewton(initialGuess, massConservation, massConservationDerivative, 1e-8);
            }
            catch (const std::exception& exc)
            {
                std::cerr << "\n\n Caught exception in computeGasPlumeDist, the local non-linear solver did not converge! Maybe the initialGuess equals the domainHeight, which is a problem for convergence?" << exc.what() << std::endl;
                throw;
            }
        }

        //check if gas plume distance is within domain
        if (gasPlumeDistance<0.0 || gasPlumeDistance > domainHeight)
        {
            std::ostringstream message;
            message << "Gas plume distance is outside the column: " << gasPlumeDistance << " not in [0, " << domainHeight << "]. " << "Wetting-phase saturation: " << satWCoarse << ", densities: " << densityW << ", " << densityNw;
            throw std::runtime_error(message.str());
        }

        return gasPlumeDistance;
    }


    /*!
     * \brief Computes the capillary pressure on the coarse level
     *
     * \param gasPlumeDist    value of the gas plume distance
     * \param densities       contains the phase densities (here: 2 phases)
     * \param gravityNorm     norm of the gravity
     * \param entryPressureBC entry pressure of the Brook-Corey model
     */
    const Scalar computeCapillaryPressureCoarse(const Scalar& gasPlumeDist,
                                                const PhaseDensities& densities,
                                                const Scalar& gravityNorm,
                                                const Scalar& entryPressureBC) const
    {
        const Scalar referenceDensityW = densities.wetting;
        const Scalar referenceDensityNw = densities.nonwetting;

        //calculate the coarse-level capillary pressure
        const Scalar pcCoarse = gravityNorm * gasPlumeDist * (referenceDensityNw - referenceDensityW) + entryPressureBC;
        return pcCoarse;
    }


    /*!
     * \brief Reconstructs the fine-level wetting-phase and non-wetting phase pressures

     *
     * \param gasPlumeDist      value of the gas plume distance
     * \param densities         contains the phase densities (here: 2 phases)
     * \param gravityNorm       norm of the gravity
     * \param heightAboveBottom height (relative to the bottom of the domain) at which the pressures should be reconstructed
     * \param pressureWCoarse   wetting-phase pressure on the coarse level
     * \param entryPressureBC   entry pressure of the Brook-Corey model
     */
    const std::vector<Scalar> reconstPressure(const Scalar& gasPlumeDist,
                                              const PhaseDensities& densities,
                                              const Scalar& gravityNorm,
                                              const Scalar& heightAboveBottom,
                                              const Scalar& pressureWCoarse,
                                              const Scalar& entryPressureBC) const
    {
        const Scalar referenceDensityW = densities.wetting;
        const Scalar referenceDensityNw = densities.nonwetting;

        std::vector<Scalar> reconstructedPressures(numPhases, 0.0);

        //So far, only the capillaryFringe model is implemented
        if(heightAboveBottom <= gasPlumeDist)
        {
            reconstructedPressures[waterPhaseIdx] = pressureWCoarse - referenceDensityW * gravityNorm * heightAboveBottom;
            reconstructedPressures[gasPhaseIdx] = reconstructedPressures[waterPhaseIdx] + entryPressureBC;
        }
        else
        {
            reconstructedPressures[waterPhaseIdx] = pressureWCoarse - referenceDensityW * gravityNorm * heightAboveBottom;
            reconstructedPressures[gasPhaseIdx] = pressureWCoarse - referenceDensityW * gravityNorm * gasPlumeDist - referenceDensityNw * gravityNorm * (heightAboveBottom - gasPlumeDist) + entryPressureBC;
        }

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
     * \param deltaZ                discretizaion width of the fine-level grid in vertical direction
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     */
    const Scalar reconstructSaturation(const GasPlumeDistances& gasPlumeDistances,
                                       const PhaseDensities& densities,
                                       const ResidualSaturations& residualSaturations,
                                       const Scalar& gravityNorm,
                                       const Scalar& heightAboveBottom,
                                       const Scalar& deltaZ,
                                       const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar gasPlumeDist = gasPlumeDistances.current;
        const Scalar minGasPlumeDist = gasPlumeDistances.minimum;

        const Scalar lowerBound = heightAboveBottom - deltaZ/2.0;
        const Scalar upperBound = heightAboveBottom + deltaZ/2.0;
        const Scalar cellHeight = upperBound - lowerBound;

        const Scalar targetAverageError = 1e-10;
        const Scalar targetIntegralError = targetAverageError * cellHeight;

        const Scalar snr = residualSaturations.nonwetting;
        const Scalar saturationWBelowPlume = 1.0;
        const Scalar saturationWTrappedRegion = 1.0 - snr;
        const auto saturationWAbovePlume = [&](Scalar z)
        {
            return saturationWAbovePlume_(z, gasPlumeDist, densities, residualSaturations, gravityNorm, brooksCoreyParameters);
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
            saturationIntegral = integrateScalarFunction(saturationWAbovePlume, lowerBound, upperBound, targetIntegralError);
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
            const Scalar integralAbovePlume = integrateScalarFunction(saturationWAbovePlume, gasPlumeDist, upperBound, targetIntegralError);
            saturationIntegral = integralTrappedRegion + integralAbovePlume;
        }
        else // 6)
        {
            const Scalar integralBelowMinimum = saturationWBelowPlume * (minGasPlumeDist - lowerBound);
            const Scalar integralTrappedRegion = saturationWTrappedRegion * (gasPlumeDist - minGasPlumeDist);
            const Scalar integralAbovePlume = integrateScalarFunction(saturationWAbovePlume, gasPlumeDist, upperBound, targetIntegralError);
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
     * \param entryPressureBC   entry pressure of the Brook-Corey model
     */
    const Scalar reconstCapillaryPressure(const Scalar& gasPlumeDist,
                                          const PhaseDensities& densities,
                                          const Scalar& gravityNorm,
                                          const Scalar& heightAboveBottom,
                                          const Scalar& entryPressureBC) const
    {
        const Scalar referenceDensityW = densities.wetting;
        const Scalar referenceDensityNw = densities.nonwetting;

        Scalar reconstCapillaryPressure = 0.0;

        //for capillary fringe model
        if(heightAboveBottom <= gasPlumeDist)
        {
            reconstCapillaryPressure = entryPressureBC;
        }
        else if(heightAboveBottom > gasPlumeDist)
        {
            reconstCapillaryPressure = referenceDensityW * gravityNorm * (heightAboveBottom - gasPlumeDist) + entryPressureBC - referenceDensityNw * gravityNorm * (heightAboveBottom - gasPlumeDist);
        }

        return reconstCapillaryPressure;
    }


    /*!
     * \brief Reconstructs the fine-level wetting-phase and non-wetting-phase mobilities
     *
     * \param gasPlumeDistances     contains the gas plume distance and the minimum gas plume distance
     * \param densities             contains the phase densities (here: 2 phases)
     * \param viscosities           contains the phase viscosities (here: 2 phases)
     * \param residualSaturations   contains the phase residual saturation (here: 2 phases)
     * \param gravityNorm           norm of the gravity
     * \param heightAboveBottom     height (relative to the bottom of the domain) at which the mobilities should be reconstructed
     * \param deltaZ                discretizaion width of the fine-level grid in vertical direction
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     */
    const std::vector<Scalar> reconstMobilitiesFine(const GasPlumeDistances& gasPlumeDistances,
                                                    const PhaseDensities& densities,
                                                    const PhaseViscosities& viscosities,
                                                    const ResidualSaturations& residualSaturations,
                                                    const Scalar& gravityNorm,
                                                    const Scalar& heightAboveBottom,
                                                    const Scalar& deltaZ,
                                                    const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar gasPlumeDist = gasPlumeDistances.current;
        const Scalar viscosityW = viscosities.wetting;
        const Scalar viscosityNw = viscosities.nonwetting;

        std::vector<Scalar> mobilitesFine(numPhases);

        const Scalar lowerBound = heightAboveBottom - deltaZ/2.0;
        const Scalar upperBound = heightAboveBottom + deltaZ/2.0;

        Scalar mobilityWFine = 0.0;
        Scalar mobilityNwFine = 0.0;

        const Scalar relPermWBelowPlume = 1.0;  //constant value
        const Scalar relPermNwBelowPlume = 0.0; //constant value

        const auto saturationWAbovePlume = [&](Scalar z)
        {
            return saturationWAbovePlume_(z, gasPlumeDist, densities, residualSaturations, gravityNorm, brooksCoreyParameters);
        };

        //define function for wetting-phase relative peremability (here Brooks-Corey)
        const auto relPermWAbovePlume = [&saturationWAbovePlume, &brooksCoreyParameters, &residualSaturations](Scalar z)
        {
            auto [swr, snr] = residualSaturations;
            auto lambdaBC = brooksCoreyParameters.lambda;
            Scalar swe = (saturationWAbovePlume(z)-swr)/(1-swr-snr);

            return std::pow(swe, 2.0/lambdaBC + 3.0);
        };

        //define function for non-wetting-phase relative peremability (here Brooks-Corey)
        const auto relPermNwAbovePlume = [&saturationWAbovePlume, &brooksCoreyParameters, &residualSaturations](Scalar z)
        {
            auto [swr, snr] = residualSaturations;
            auto lambdaBC = brooksCoreyParameters.lambda;

            Scalar swe = (saturationWAbovePlume(z)-swr)/(1-swr-snr);
            const Scalar exponent = 2.0/lambdaBC + 1.0;
            const Scalar sne = 1.0 - swe;
            return sne*sne*(1.0 - std::pow(swe, exponent));
        };

        // splitting the integration interval explicitly at gasPlumeDist is cheaper for the integrator, since integrand is piecewise at gasPlumeDist. If snr!=0, there is also a jump in the saturation at gasPlumeDist
        const auto integrateRelativePermeability = [&](const auto& integrandAbovePlume, const Scalar relPermBelowPlume)
        {
            const Scalar cellHeight = upperBound - lowerBound;
            const Scalar targetAverageError = 1e-10;
            const Scalar targetIntegralError = targetAverageError * cellHeight; // since the integral is divided by cellHeight later on for the computation of the mobility

            // differentiate between three cases:
            //     1) fine-level element is completely under the gas plume,
            //     2) fine-level element is completely above the gas plume,
            //     3) fine-level element intersects with gas plume
            //     for zp_min<z<zp, there might be snr present and the rel perm for gas is treated as 0, thus the check z<zp is enough
            if (upperBound <= gasPlumeDist)
                return relPermBelowPlume * cellHeight;
            else if (lowerBound >= gasPlumeDist)
                return integrateScalarFunction(integrandAbovePlume, lowerBound, upperBound, targetIntegralError);
            else
            {
                // part of cell below gas plume
                const Scalar integralBelowPlume = relPermBelowPlume * (gasPlumeDist - lowerBound);

                // part of cell above gas plume
                const Scalar integralAbovePlume = integrateScalarFunction(integrandAbovePlume, gasPlumeDist, upperBound, targetIntegralError);

                return integralBelowPlume  + integralAbovePlume;
            }
        };

        mobilityWFine = integrateRelativePermeability(relPermWAbovePlume, relPermWBelowPlume);
        mobilityNwFine = integrateRelativePermeability(relPermNwAbovePlume, relPermNwBelowPlume);

        //average
        mobilityWFine = mobilityWFine/(upperBound - lowerBound);
        mobilityNwFine = mobilityNwFine/(upperBound - lowerBound);

        //turn relative permeability to mobility
        mobilityWFine = mobilityWFine/viscosityW;
        mobilityNwFine = mobilityNwFine/viscosityNw;

        //store
        mobilitesFine[waterPhaseIdx] = mobilityWFine;
        mobilitesFine[gasPhaseIdx] = mobilityNwFine;

        return mobilitesFine;
    }

private:

    /*!
     * \brief Helper function for evaluation wetting-phase saturation above the gas plume
     *
     * \param heightAboveBottom     height (relative to the bottom of the domain) at which the saturation should be reconstructed
     * \param gasPlumeDist          gas plume distance
     * \param densities             contains the phase densities (here: 2 phases)
     * \param residualSaturations   contains the phase residual saturation (here: 2 phases)
     * \param gravityNorm           norm of the gravity
     * \param brooksCoreyParameters contains the two Brooks-Corey parameters (lambda and entry pressure)
     */
    Scalar saturationWAbovePlume_(Scalar heightAboveBottom,
                                  Scalar gasPlumeDist,
                                  const PhaseDensities& densities,
                                  const ResidualSaturations& residualSaturations,
                                  Scalar gravityNorm,
                                  const BrooksCoreyParameters& brooksCoreyParameters) const
    {
        const Scalar referenceDensityW = densities.wetting;
        const Scalar referenceDensityNw = densities.nonwetting;
        const Scalar swr = residualSaturations.wetting;
        const Scalar snr = residualSaturations.nonwetting;
        const Scalar lambdaBC = brooksCoreyParameters.lambda;
        const Scalar entryPressureBC = brooksCoreyParameters.entryPressure;

        const Scalar satW = std::pow((entryPressureBC + (referenceDensityW-referenceDensityNw)*gravityNorm*(heightAboveBottom-gasPlumeDist)),-lambdaBC)*std::pow(entryPressureBC, lambdaBC)*(1-swr-snr)+swr;

        return satW;
    }

    bool float_equal_(Scalar a,
                      Scalar b,
                      Scalar epsilon = std::numeric_limits<Scalar>::epsilon()) const
    {
        return std::abs(a - b) <= epsilon * std::max(std::abs(a), std::abs(b));
    }
};

} // end namespace Dumux

#endif
