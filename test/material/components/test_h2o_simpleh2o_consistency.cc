// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief Test some components for consistency of their physical properties.
 */

#include "config.h"

#include <cmath>
#include <dune/common/exceptions.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/simpleh2o.hh>

int main(int argc, char *argv[])
{
    using namespace Dumux;

    // See the currently allowed range in simpleh2O.hh
    constexpr double lowerTemperatureBound = 273.15 + 0.0; // 0°C
    constexpr double upperTemperatureBound = 273.15 + 150.0; // 150°C
    constexpr int nTemperatures = 20;
    constexpr double pressure = 1.0e5; // 1 bar

    for (int i = 0; i < nTemperatures; ++i)
    {
        const double temperature = lowerTemperatureBound + (upperTemperatureBound - lowerTemperatureBound) * i / (nTemperatures - 1);
        const double simpleH2OGasEnthalpy = Components::SimpleH2O<double>::gasEnthalpy(temperature, pressure);
        const double h2OGasEnthalpy = Components::H2O<double>::gasEnthalpy(temperature, pressure);
        // compare the two gas enthalpies
        const double relDiff = (simpleH2OGasEnthalpy - h2OGasEnthalpy) / h2OGasEnthalpy;
        const double relTol = 0.05; // We allow a relative error of 5%
        if (std::abs(relDiff) > relTol) {
            DUNE_THROW(Dune::Exception, "The gas enthalpy of the simple H2O component is not consistent with the H2O component. "
                                        << "The relative difference " << relDiff << " at temperature " << temperature << " K."
                                        << " is larger than the relative tolerance " << relTol << ".");
        }
    }
    return 0;
}
