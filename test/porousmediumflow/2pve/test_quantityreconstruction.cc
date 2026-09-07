// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief Unit test for quantity reconstructor.
 */

#include <config.h>
#include <initializer_list>
#include <iomanip>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

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

    template<class F>
    void expectThrow(F&& function, const std::string& description)
    {
        try
        {
            function();
        }
        catch (const Dune::InvalidStateException&)
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
                               brooksCoreyParameters);
    constexpr int numCells = 100;
    const Scalar cellHeight = domainHeight/numCells;
    Scalar reconstructedSwAverage = 0.0;
    Scalar mobilityWCoarse = 0.0;
    Scalar mobilityNwCoarse = 0.0;
    const Scalar permeabilityFine = 2e-12;
    const Scalar permeabilityCoarse = permeabilityFine*domainHeight;
    for (int cellIdx = 0; cellIdx < numCells; ++cellIdx)
    {
        const Scalar cellCenter = (cellIdx + 0.5)*cellHeight;
        // assumes that porosity is constant everywhere
        const Scalar reconstructedSw = reconstructor.reconstructSaturation(
                                        GasPlumeDistances{computedZp, computedZp},
                                        densities,
                                        residualSaturations,
                                        gravity,
                                        cellCenter,
                                        cellHeight,
                                        brooksCoreyParameters);
        reconstructedSwAverage += reconstructedSw*cellHeight;

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
    reconstructedSwAverage /= domainHeight;
    TwoPVE::checkClose(reconstructedSwAverage, coarseSaturationW, 1.0e-8, "coarse-level wetting-phase saturation");

    mobilityWCoarse /= permeabilityCoarse;
    mobilityNwCoarse /= permeabilityCoarse;
    // regression values for the column-integrated mobilities. These values protect against unintended changes in reconstruction and upscaling
    const Scalar mobilityWRef = 498.59507631;
    const Scalar mobilityNwRef = 7710.08376616;
    TwoPVE::checkClose(mobilityWCoarse, mobilityWRef, 1.0e-8, "coarse-level wetting-phase mobility");
    TwoPVE::checkClose(mobilityNwCoarse, mobilityNwRef, 1.0e-8, "coarse-level non-wetting-phase mobility");

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
            brooksCoreyParameters
        );
    };
    TwoPVE::expectThrow(zpCallGravity, "zero gravity");

    auto zpCallEqualDensities = [&]()
    {
        reconstructor.computeGasPlumeDist(
            PhaseDensities{1000.0, 1000.0},
            residualSaturations,
            gravity,
            domainHeight,
            coarseSaturationW,
            brooksCoreyParameters
        );
    };
    TwoPVE::expectThrow(zpCallEqualDensities, "equal phase densities");

    auto zpCallDensities = [&]()
    {
        reconstructor.computeGasPlumeDist(
            PhaseDensities{900.0, 1000.0},
            residualSaturations,
            gravity,
            domainHeight,
            coarseSaturationW,
            brooksCoreyParameters
        );
    };
    TwoPVE::expectThrow(zpCallDensities, "wetting-phase density that is smaller than the non-wetting one.");

    auto zpCallLambda = [&]()
    {
        reconstructor.computeGasPlumeDist(
            densities,
            residualSaturations,
            gravity,
            domainHeight,
            coarseSaturationW,
            BrooksCoreyParameters{1.0, 1.0e5}
        );
    };
    TwoPVE::expectThrow(zpCallLambda, "Brooks-Corey lambda equal to one");

    return 0;
}
