//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <vector>

#include <dune/istl/bvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/initialize.hh>
#include <dumux/nonlinear/findscalarroot.hh>
#include <dumux/io/gnuplotinterface.hh>

#include <dumux/experimental/timestepping/newmarkbeta.hh>

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // integrate d^2x/dt^2 = -x with Newmark-beta scheme
    Experimental::NewmarkBeta<double, Dune::BlockVector<double>> newmark;
    Dune::BlockVector<double> x(1); x = -1.0;
    Dune::BlockVector<double> v(1); v = 0.0;
    Dune::BlockVector<double> a(1); a = 1.0;
    newmark.initialize(x, v, a);

    const double dt = 0.1;
    const double tEnd = 16*M_PI;
    const int nSteps = tEnd/dt;

    const bool showPlot = getParam<bool>("Plot", false);
    Dumux::GnuplotInterface<double> plotter(showPlot);
    std::vector<double> tValues, xValues, xExact, vValues, vExact, aValues, aExact;
    tValues.reserve(nSteps);
    xValues.reserve(nSteps); xExact.reserve(nSteps);
    vValues.reserve(nSteps); vExact.reserve(nSteps);
    aValues.reserve(nSteps); aExact.reserve(nSteps);

    for (int i = 0; i < nSteps; ++i)
    {
        const double t = dt*i;

        // mimic a nonlinear solver by finding the root of the residual (d^2x/dt^2 + x = 0)
        const auto residual = [&](const auto x){ return newmark.acceleration(0, dt, x) + x; };
        const auto xNew = findScalarRootNewton(x[0], residual);

        x[0] = xNew;
        // update with the new position
        // after the update we have the current position, velocity and acceleration
        newmark.update(dt, x);

        tValues.push_back(t);
        xValues.push_back(newmark.position(0));
        xExact.push_back(-std::cos(t));
        vValues.push_back(newmark.velocity(0));
        vExact.push_back(std::sin(t));
        aValues.push_back(newmark.acceleration(0));
        aExact.push_back(std::cos(t));

        if (i % 50 == 0)
        {
            plotter.resetPlot();
            plotter.addDataSetToPlot(tValues, xValues, "x.dat", "w l");
            plotter.addDataSetToPlot(tValues, xExact, "xe.dat", "w l");
            plotter.addDataSetToPlot(tValues, vValues, "v.dat", "w l");
            plotter.addDataSetToPlot(tValues, vExact, "ve.dat", "w l");
            plotter.addDataSetToPlot(tValues, aValues, "a.dat", "w l");
            plotter.addDataSetToPlot(tValues, aExact, "ae.dat", "w l");
            plotter.plot("result");
        }
    }

    // compare the final position with the exact solution
    const double error = std::abs(x[0] + std::cos(tEnd));
    if (error > 0.006)
        DUNE_THROW(Dune::Exception, "Integration error too large: " << error);

    return 0;
}
