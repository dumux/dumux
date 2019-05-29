// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Test the spline law
 */
#include <config.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscoreyparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/splinelaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/splinelawparams.hh>

std::vector<double> linspace(const double begin, const double end, const double samples)
{
    const double delta = (end-begin)/static_cast<double>(samples-1);
    std::vector<double> vec(samples);
    for (int i = 0; i < samples; ++i)
        vec[i] = begin + i*delta;
    return vec;
}

template<class Function>
std::vector<double> eval(const Function& f, const std::vector<double>& x)
{
    auto y = x;
    std::transform(x.begin(), x.end(), y.begin(), [&](const double x) { return f(x); });
    return y;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    using Scalar = double;
    using namespace Dumux;
    RegularizedBrooksCoreyParams<Scalar> bcParams(1e4, 2.0);

    // we test the spline interpolation against a sample function
    const auto pc = [&](const double x){ return RegularizedBrooksCorey<Scalar>::pc(bcParams, x); };
    const auto dpc = [&](const double x){ return RegularizedBrooksCorey<Scalar>::dpc_dswe(bcParams, x); };

    // create some test samples
    const auto testPoints = linspace(0.0, 1.0, 1000);
    const auto ref = eval(pc, testPoints);
    const auto refDeriv = eval(dpc, testPoints);

    // create the spline sample points
    const std::vector<double> swData = {0.0, 0.05, 0.2, 0.7, 1.0};
    const auto pcData = eval(pc, swData);

    SplineLawParams<Scalar> slParams(swData, pcData, pcData, pcData);

    // evaluate spline law and derivative
    const auto result = eval([&](const double x) { return SplineLaw<Scalar>::pc(slParams, x); }, testPoints);
    const auto resultDeriv = eval([&](const double x) { return SplineLaw<Scalar>::dpc_dswe(slParams, x); }, testPoints);

    // plot with Gnuplot
    Dumux::GnuplotInterface<double> gnuplot(/*persist=*/true);
    gnuplot.addDataSetToPlot(testPoints, ref, "pc reference");
    gnuplot.addDataSetToPlot(testPoints, result, "pc spline law");
    // gnuplot.addDataSetToPlot(testPoints, refDeriv, "dpc reference");
    // gnuplot.addDataSetToPlot(testPoints, resultDeriv, "dpc spline law");
    gnuplot.plot("spline law");

    return 0;
}
