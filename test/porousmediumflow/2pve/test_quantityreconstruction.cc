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

#include <cmath>

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
        if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::CmpStyle::absolute>(actual, expected, tolerance))
            DUNE_THROW(Dune::Exception, "Unexpected " << quantity << ": expected "
                       << expected << ", obtained " << actual);
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

    // test coarse-level saturation
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
    }
    reconstructedSwAverage /= domainHeight;
    TwoPVE::checkClose(reconstructedSwAverage, coarseSaturationW, 1.0e-8, "column-averaged reconstructed saturation");

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
