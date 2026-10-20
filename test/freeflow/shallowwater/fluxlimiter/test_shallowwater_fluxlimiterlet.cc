// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the mobility of the LET flux limiter for small water depths:
 *        zero at and below the lower depth, one at and above the upper depth,
 *        continuous and increasing in between
 */
#include <config.h>

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/flux/shallowwater/fluxlimiterlet.hh>

int main()
{
    using namespace Dumux;

    for (const auto [lowerH, upperH] : {std::pair{1e-5, 1e-3}, std::pair{0.0, 1e-2}, std::pair{2e-3, 5e-3}})
    {
        auto mobility = [&](double h) { return ShallowWater::fluxLimiterLET(h, h, upperH, lowerH); };

        if (mobility(lowerH) != 0.0 || mobility(0.5*lowerH) != 0.0 || mobility(0.0) != 0.0)
            DUNE_THROW(Dune::Exception, "Mobility is not zero at and below the lower water depth " << lowerH);
        if (mobility(upperH) != 1.0 || mobility(2.0*upperH) != 1.0)
            DUNE_THROW(Dune::Exception, "Mobility is not one at and above the upper water depth " << upperH);
        if (std::abs(mobility(0.5*(lowerH + upperH)) - 0.5) > 1e-14)
            DUNE_THROW(Dune::Exception, "Mobility is not one half at the mean of the lower and upper water depth");

        constexpr int n = 1000;
        double previous = 0.0;
        for (int i = 1; i <= n; ++i)
        {
            const double h = lowerH + (upperH - lowerH)*i/n;
            const double m = mobility(h);
            if (!(m >= previous) || m - previous > 0.01)
                DUNE_THROW(Dune::Exception, "Mobility is not continuous and increasing at h = " << h);
            previous = m;
        }
        std::cout << "LET mobility checks passed for lower depth " << lowerH << " and upper depth " << upperH << std::endl;
    }

    return 0;
}
