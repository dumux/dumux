// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Test for the analytical solution of Mandel's problem
 *
 * The series solution at \f$t = 0\f$ has to reproduce the undrained response,
 * the pressure vanishes at the drained boundary and exceeds the undrained
 * pressure at the center for early times (Mandel-Cryer effect).
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include "analyticalsolution.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using GlobalPosition = Dune::FieldVector<double, 2>;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    const MandelAnalyticalSolution<double> solution;
    const auto upperRight = getParam<GlobalPosition>("Grid.UpperRight");
    const double a = upperRight[0];
    const double b = upperRight[1];
    const double pInit = solution.initialPressure(upperRight);
    const auto uInit = solution.initialDisplacement(upperRight);

    using std::abs;
    for (const double xRel : {0.0, 0.25, 0.5, 0.75})
    {
        const GlobalPosition pos{xRel*a, b};

        const auto u = solution.displacement(pos, 0.0);
        const auto u0 = solution.initialDisplacement(pos);
        for (int dir = 0; dir < 2; ++dir)
            if (abs(u[dir] - u0[dir]) > 1e-3*abs(uInit[dir]))
                DUNE_THROW(Dune::Exception, "Displacement " << dir << " at t = 0, x/a = " << xRel
                           << ": " << u[dir] << " != " << u0[dir]);

        const double p = solution.pressure(pos, 0.0);
        if (abs(p - pInit) > 1e-2*pInit)
            DUNE_THROW(Dune::Exception, "Pressure at t = 0, x/a = " << xRel << ": " << p << " != " << pInit);
    }

    for (const double t : {0.0, 100.0, 1000.0})
        if (abs(solution.pressure(upperRight, t)) > 1e-12*pInit)
            DUNE_THROW(Dune::Exception, "Pressure at the drained boundary at t = " << t << ": " << solution.pressure(upperRight, t));

    const GlobalPosition center{0.0, 0.0};
    if (!(solution.pressure(center, 1000.0) > pInit))
        DUNE_THROW(Dune::Exception, "No Mandel-Cryer effect: p(0, 1000 s) = " << solution.pressure(center, 1000.0)
                   << " <= " << pInit);

    std::cout << "All checks passed" << std::endl;
    return 0;
}
