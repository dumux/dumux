// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Test for the outer cell states of the shallow water boundary conditions
 */
#include <config.h>

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>

namespace {

using Scalar = double;
using Normal = Dune::FieldVector<Scalar, 2>;

void check(bool condition, const std::string& what)
{
    if (!condition)
        DUNE_THROW(Dune::Exception, "Check failed: " << what);
}

bool close(Scalar a, Scalar b, Scalar tol = 1e-12)
{ return Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(a, b, tol); }

} // end namespace

int main(int argc, char** argv)
{
    using namespace Dumux;

    const Scalar gravity = 9.81;
    const Normal outward{1.0, 0.0};
    const Normal inward{-1.0, 0.0};
    const Normal side{0.0, 1.0};

    // slip wall: normal velocity mirrored, tangential velocity kept, no mass flux
    {
        const auto state = ShallowWater::wallBoundary(0.2, 1.5, 0.3, side);
        check(close(state[0], 0.2) && close(state[1], 1.5) && close(state[2], -0.3),
              "wall mirrors the normal velocity only");

        const auto flux = ShallowWater::riemannProblem(0.2, state[0], 1.5, state[1], 0.3, state[2],
                                                       0.0, 0.0, gravity, side);
        check(close(flux[0], 0.0), "wall carries no mass flux");
    }

    // inflow: prescribed depth and speed, directed into the domain
    {
        const auto state = ShallowWater::inflowBoundary(0.05, 2.0, inward);
        check(close(state[0], 0.05) && close(state[1], 2.0) && close(state[2], 0.0),
              "inflow state points into the domain");
    }

    // free overfall: critical depth for subcritical interior, pass-through for supercritical
    {
        const Scalar h = 0.1, u = 0.5;
        const auto state = ShallowWater::criticalDepthOutflowBoundary(h, u, 0.0, gravity, outward);
        const Scalar criticalDepth = std::cbrt(h*u*h*u/gravity);
        check(close(state[0], criticalDepth), "subcritical outflow reaches critical depth");
        check(close(state[1]*state[1], gravity*criticalDepth) && state[1] > 0.0,
              "critical outflow moves at the wave celerity");
        check(close(state[0]*state[1], h*u), "critical outflow keeps the unit discharge");

        const auto supercritical = ShallowWater::criticalDepthOutflowBoundary(0.05, 3.0, 0.0, gravity, outward);
        check(close(supercritical[0], 0.05) && close(supercritical[1], 3.0),
              "supercritical outflow copies the interior state");

        const Scalar v = 0.3;
        const auto oblique = ShallowWater::criticalDepthOutflowBoundary(h, u, v, gravity, outward);
        check(close(oblique[0], criticalDepth) && close(oblique[2], v),
              "oblique outflow reaches critical depth for the normal discharge and keeps the tangential velocity");

        // the speed is supercritical, but the normal velocity is not
        const Scalar hFast = 0.05, uNormal = 0.3, uTangential = 3.0;
        const auto tangential = ShallowWater::criticalDepthOutflowBoundary(hFast, uNormal, uTangential, gravity, outward);
        check(close(tangential[0], std::cbrt(hFast*uNormal*hFast*uNormal/gravity)) && close(tangential[2], uTangential),
              "the flow regime at the boundary is determined by the normal velocity");

        const auto backflow = ShallowWater::criticalDepthOutflowBoundary(h, -u, v, gravity, outward);
        const auto wall = ShallowWater::wallBoundary(h, -u, v, outward);
        check(close(backflow[0], wall[0]) && close(backflow[1], wall[1]) && close(backflow[2], wall[2]),
              "water flowing towards the domain at a free overfall is stopped by a wall");
    }

    // vanishing discharge behaves as a wall
    {
        const auto wall = ShallowWater::wallBoundary(0.2, 0.4, 0.1, inward);
        const auto state = ShallowWater::fixedDischargeBoundary(0.0, 0.2, 0.4, 0.1, gravity, inward);
        check(close(state[0], wall[0]) && close(state[1], wall[1]) && close(state[2], wall[2]),
              "zero discharge gives the wall state");
    }

    // fixed discharge: the outer state carries the prescribed discharge through the face
    {
        const Scalar discharge = -0.05;
        const auto state = ShallowWater::fixedDischargeBoundary(discharge, 0.2, 0.2, 0.0, gravity, inward);
        const auto normalDischarge = state[0]*(inward[0]*state[1] + inward[1]*state[2]);
        check(close(normalDischarge, discharge, 1e-9), "outer state carries the prescribed discharge");
    }

    std::cout << "All shallow water boundary state checks passed" << std::endl;
    return 0;
}
