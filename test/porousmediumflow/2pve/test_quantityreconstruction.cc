// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Unit test for quantity reconstructor.
 */

#include <config.h>
#include <algorithm>
#include <initializer_list>
#include <iomanip>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/integrate.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/porousmediumflow/2pve/quantityreconstruction.hh>

#include "properties.hh"

namespace TwoPVE {

    template<class Scalar>
    void checkClose(const Scalar actual,
                    const Scalar expected,
                    const Scalar tolerance,
                    const std::string& quantity)
    {
        int digits = static_cast<int>(std::ceil(-std::log10(tolerance)));
        if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::CmpStyle::absolute>(actual, expected, tolerance))
            DUNE_THROW(Dune::Exception, std::fixed << std::setprecision(digits) << "Unexpected " << quantity << ": expected " << expected << ", obtained " << actual);
    }

    template<class Exception, class F>
    void expectThrow(F&& function, const std::string& description)
    {
        try
        {
            function();
        }
        catch (const Exception&)
        {
            return;
        }

        DUNE_THROW(Dune::Exception, "Expected an exception for " << description);
    }
} // end namespace TwoPVE

int main()
{
    using TypeTag = Dumux::Properties::TTag::TwoPVEImmiscibleTpfa;
    using Scalar = Dumux::GetPropType<TypeTag, Dumux::Properties::Scalar>;
    using Reconstructor = Dumux::TwoPVEQuantityReconst<TypeTag>;
    using GasPlumeDistances = Dumux::TwoPVE::GasPlumeDistancesData<Scalar>;
    using PhaseDensities = Dumux::TwoPVE::PhaseDensitiesData<Scalar>;
    using PhaseViscosities = Dumux::TwoPVE::PhaseViscositiesData<Scalar>;
    using ResidualSaturations = Dumux::TwoPVE::ResidualSaturationsData<Scalar>;
    using BrooksCoreyParameters = Dumux::TwoPVE::BrooksCoreyParametersData<Scalar>;

    const Reconstructor reconstructor;
    const PhaseDensities densities{1000.0, 100.0};
    const PhaseViscosities viscosities{1.0e-3, 1.0e-5};
    const ResidualSaturations residualSaturations{0.1, 0.2};
    const BrooksCoreyParameters brooksCoreyParameters{2.0, 1.0e5};
    const Scalar gravity = 9.81;
    const Scalar domainHeight = 10.0;

    // test gas plume distance
    const auto gasPlumeDistanceNoGas = reconstructor.computeGasPlumeDist(
                                        densities,
                                        residualSaturations,
                                        gravity,
                                        domainHeight,
                                        1.0,
                                        domainHeight,
                                        brooksCoreyParameters);
    TwoPVE::checkClose(gasPlumeDistanceNoGas, domainHeight, 1.0e-12, "gas plume distance for a fully water-saturated column");

    // test reconstructed capillary pressure
    const Scalar gasPlumeDistance = 4.0;
    const Scalar pressureWCoarse = 2.0e5;
    for (const Scalar height : {2.0, gasPlumeDistance, 7.0})
    {
        const auto pressures = reconstructor.reconstPressure(
                                gasPlumeDistance,
                                densities,
                                gravity,
                                height,
                                pressureWCoarse,
                                brooksCoreyParameters.entryPressure);
        const auto capillaryPressure = reconstructor.reconstCapillaryPressure(
                                        gasPlumeDistance,
                                        densities,
                                        gravity,
                                        height,
                                        brooksCoreyParameters.entryPressure);
        TwoPVE::checkClose(pressures[1] - pressures[0], capillaryPressure, 1.0e-10, "pressure/capillary-pressure consistency");
    }

    // test reconstructed wetting-phase saturations
    // zp equals short notation of gas plume distance
    const auto zpMin = 2.0;
    const auto heightBelowZp = 1.0;
    const auto deltaZ = 1.0;
    const auto reconstructedSwBelowZpMin = reconstructor.reconstructSaturation(
                                            GasPlumeDistances{gasPlumeDistance, zpMin},
                                            densities,
                                            residualSaturations,
                                            gravity,
                                            heightBelowZp,
                                            deltaZ,
                                            brooksCoreyParameters);
    TwoPVE::checkClose(reconstructedSwBelowZpMin, 1.0, 1.0e-12, "saturation below the historical plume");
    const auto heightTrappedGasRegion = 3.0;
    const auto reconstructedSwTrappedGasRegion = reconstructor.reconstructSaturation(
                                                  GasPlumeDistances{gasPlumeDistance, zpMin},
                                                  densities,
                                                  residualSaturations,
                                                  gravity,
                                                  heightTrappedGasRegion,
                                                  deltaZ,
                                                  brooksCoreyParameters);
    TwoPVE::checkClose(reconstructedSwTrappedGasRegion, 1.0 - residualSaturations.nonwetting, 1.0e-12,
               "saturation in the trapped-gas region");

    // test reconstructed mobilities
    const auto mobilities = reconstructor.reconstMobilitiesFine(
                             GasPlumeDistances{gasPlumeDistance, gasPlumeDistance},
                             densities,
                             viscosities,
                             residualSaturations,
                             gravity,
                             heightBelowZp,
                             deltaZ,
                             brooksCoreyParameters);
    const auto relPermWBelowZp = 1.0;
    TwoPVE::checkClose(mobilities[0], relPermWBelowZp/viscosities.wetting, 1.0e-10, "wetting mobility below the plume");
    TwoPVE::checkClose(mobilities[1], 0.0, 1.0e-12, "nonwetting mobility below the plume");

    // test coarse-level saturation and mobility
    const Scalar coarseSaturationW = 0.7;
    const Scalar computedZp = reconstructor.computeGasPlumeDist(
                               densities,
                               residualSaturations,
                               gravity,
                               domainHeight,
                               coarseSaturationW,
                               domainHeight,
                               brooksCoreyParameters);
    constexpr int numCells = 100;
    const Scalar cellHeight = domainHeight/numCells;

    // assumes that porosity is constant everywhere
    const auto columnAverageSaturationW = [&](const Scalar zp, const Scalar zpMin, const BrooksCoreyParameters& materialParameters)
    {
        Scalar saturationIntegral = 0.0;
        for (int cellIdx = 0; cellIdx < numCells; ++cellIdx)
            saturationIntegral += reconstructor.reconstructSaturation(
                                      GasPlumeDistances{zp, zpMin},
                                      densities,
                                      residualSaturations,
                                      gravity,
                                      (cellIdx + 0.5)*cellHeight,
                                      cellHeight,
                                      materialParameters)*cellHeight;
        return saturationIntegral/domainHeight;
    };
    TwoPVE::checkClose(columnAverageSaturationW(computedZp, computedZp, brooksCoreyParameters), coarseSaturationW, 1.0e-8, "coarse-level wetting-phase saturation");

    Scalar mobilityWCoarse = 0.0;
    Scalar mobilityNwCoarse = 0.0;
    const Scalar permeabilityFine = 2e-12;
    const Scalar permeabilityCoarse = permeabilityFine*domainHeight;
    for (int cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellCenter = (cellIdx + 0.5)*cellHeight;

        // assuming homogeneous fine-level permeability
        const std::vector<Scalar> reconstructedMobilites = reconstructor.reconstMobilitiesFine(
                                                         GasPlumeDistances{computedZp, computedZp},
                                                         densities,
                                                         viscosities,
                                                         residualSaturations,
                                                         gravity,
                                                         cellCenter,
                                                         cellHeight,
                                                         brooksCoreyParameters);
        mobilityWCoarse += permeabilityFine*reconstructedMobilites[0]*cellHeight;
        mobilityNwCoarse += permeabilityFine*reconstructedMobilites[1]*cellHeight;
    }

    mobilityWCoarse /= permeabilityCoarse;
    mobilityNwCoarse /= permeabilityCoarse;
    // regression values for the column-integrated mobilities. These values protect against unintended changes in reconstruction and upscaling
    const Scalar mobilityWRef = 498.59507631;
    const Scalar mobilityNwRef = 7710.08376616;
    TwoPVE::checkClose(mobilityWCoarse, mobilityWRef, 1.0e-8, "coarse-level wetting-phase mobility");
    TwoPVE::checkClose(mobilityNwCoarse, mobilityNwRef, 1.0e-8, "coarse-level non-wetting-phase mobility");

    // test coarse-level saturation for a column with trapped gas between the minimum gas plume distance and the gas plume distance
    const Scalar minimumZp = 2.0;
    const Scalar computedZpTrapped = reconstructor.computeGasPlumeDist(
                                      densities,
                                      residualSaturations,
                                      gravity,
                                      domainHeight,
                                      coarseSaturationW,
                                      minimumZp,
                                      brooksCoreyParameters);
    if (!(computedZpTrapped > minimumZp))
        DUNE_THROW(Dune::Exception, "Expected a trapped-gas region, obtained gas plume distance " << computedZpTrapped << " below the minimum " << minimumZp);
    TwoPVE::checkClose(columnAverageSaturationW(computedZpTrapped, minimumZp, brooksCoreyParameters), coarseSaturationW, 1.0e-8, "coarse-level wetting-phase saturation with trapped gas");

    const auto gasPlumeDistanceOnlyTrappedGas = reconstructor.computeGasPlumeDist(
                                                 densities,
                                                 residualSaturations,
                                                 gravity,
                                                 domainHeight,
                                                 1.0,
                                                 minimumZp,
                                                 brooksCoreyParameters);
    TwoPVE::checkClose(gasPlumeDistanceOnlyTrappedGas, domainHeight, 1.0e-12, "gas plume distance for a column without mobile gas");

    // test that a coarse-level saturation below the residual saturation is reported as a recoverable numerical problem
    auto zpCallInfeasibleSaturation = [&]()
    {
        reconstructor.computeGasPlumeDist(
            densities,
            residualSaturations,
            gravity,
            domainHeight,
            0.5*residualSaturations.wetting,
            domainHeight,
            brooksCoreyParameters
        );
    };
    TwoPVE::expectThrow<Dumux::NumericalProblem>(zpCallInfeasibleSaturation, "coarse-level saturation below the residual saturation");

    // test parameter constraints for gas plume distance
    const Scalar gravityZero = 0.0;
    auto zpCallGravity = [&]()
    {
        reconstructor.computeGasPlumeDist(
            densities,
            residualSaturations,
            gravityZero,
            domainHeight,
            coarseSaturationW,
            domainHeight,
            brooksCoreyParameters
        );
    };
    TwoPVE::expectThrow<Dune::InvalidStateException>(zpCallGravity, "zero gravity");

    auto zpCallEqualDensities = [&]()
    {
        reconstructor.computeGasPlumeDist(
            PhaseDensities{1000.0, 1000.0},
            residualSaturations,
            gravity,
            domainHeight,
            coarseSaturationW,
            domainHeight,
            brooksCoreyParameters
        );
    };
    TwoPVE::expectThrow<Dune::InvalidStateException>(zpCallEqualDensities, "equal phase densities");

    auto zpCallDensities = [&]()
    {
        reconstructor.computeGasPlumeDist(
            PhaseDensities{900.0, 1000.0},
            residualSaturations,
            gravity,
            domainHeight,
            coarseSaturationW,
            domainHeight,
            brooksCoreyParameters
        );
    };
    TwoPVE::expectThrow<Dune::InvalidStateException>(zpCallDensities, "wetting-phase density that is smaller than the non-wetting one.");

    // test Brooks-Corey parameters for which some integrals of powers of the effective saturation are logarithms
    for (const Scalar lambda : {1.0, 0.5})
    {
        const BrooksCoreyParameters materialParameters{lambda, brooksCoreyParameters.entryPressure};
        const Scalar zp = reconstructor.computeGasPlumeDist(
                              densities,
                              residualSaturations,
                              gravity,
                              domainHeight,
                              coarseSaturationW,
                              domainHeight,
                              materialParameters);
        TwoPVE::checkClose(columnAverageSaturationW(zp, zp, materialParameters), coarseSaturationW, 1.0e-8,
                           "coarse-level wetting-phase saturation for lambda " + std::to_string(lambda));
    }

    // test the closed-form cell averages against numerically integrated Brooks-Corey laws
    using MaterialLaw = Dumux::FluidMatrix::BrooksCoreyNoReg<Scalar>;
    const Scalar zpReference = 4.0;
    const Scalar densityDifference = densities.wetting - densities.nonwetting;
    const auto cellAverage = [&](const auto& valueAbovePlume, const Scalar valueBelowPlume, const Scalar lower, const Scalar upper)
    {
        const Scalar integralBelowPlume = valueBelowPlume*std::max(std::min(upper, zpReference) - lower, 0.0);
        const Scalar integralAbovePlume = upper > zpReference ? Dumux::integrateScalarFunction(valueAbovePlume, std::max(lower, zpReference), upper, 1.0e-14) : 0.0;
        return (integralBelowPlume + integralAbovePlume)/(upper - lower);
    };
    for (const Scalar lambda : {2.0, 1.0, 0.5})
    {
        const BrooksCoreyParameters materialParameters{lambda, brooksCoreyParameters.entryPressure};
        const MaterialLaw materialLaw(typename MaterialLaw::BasicParams(materialParameters.entryPressure, lambda),
                                     typename MaterialLaw::EffToAbsParams(residualSaturations.wetting, residualSaturations.nonwetting));
        const auto saturationW = [&](const Scalar z){ return materialLaw.sw(materialParameters.entryPressure + densityDifference*gravity*(z - zpReference)); };
        const auto relPermW = [&](const Scalar z){ return materialLaw.krw(saturationW(z)); };
        const auto relPermNw = [&](const Scalar z){ return materialLaw.krn(saturationW(z)); };

        // cells below, across and above the gas plume distance
        for (const Scalar cellCenter : {3.5, 4.25, 6.0, 9.5})
        {
            const Scalar lower = cellCenter - 0.5*deltaZ;
            const Scalar upper = cellCenter + 0.5*deltaZ;
            const std::string cell = " (lambda " + std::to_string(lambda) + ", cell center " + std::to_string(cellCenter) + ")";

            const Scalar reconstructedSw = reconstructor.reconstructSaturation(
                                               GasPlumeDistances{zpReference, zpReference},
                                               densities,
                                               residualSaturations,
                                               gravity,
                                               cellCenter,
                                               deltaZ,
                                               materialParameters);
            TwoPVE::checkClose(reconstructedSw, cellAverage(saturationW, 1.0, lower, upper), 1.0e-12, "cell-averaged saturation" + cell);

            const auto reconstructedMobilities = reconstructor.reconstMobilitiesFine(
                                                     GasPlumeDistances{zpReference, zpReference},
                                                     densities,
                                                     viscosities,
                                                     residualSaturations,
                                                     gravity,
                                                     cellCenter,
                                                     deltaZ,
                                                     materialParameters);
            TwoPVE::checkClose(reconstructedMobilities[0]*viscosities.wetting, cellAverage(relPermW, 1.0, lower, upper), 1.0e-12, "cell-averaged wetting relative permeability" + cell);
            TwoPVE::checkClose(reconstructedMobilities[1]*viscosities.nonwetting, cellAverage(relPermNw, 0.0, lower, upper), 1.0e-12, "cell-averaged nonwetting relative permeability" + cell);
        }
    }

    return 0;
}
