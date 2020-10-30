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
 * \brief Test the simple cubic spline implementation
 */
#include <config.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dumux/common/math.hh>
#include <dumux/common/monotonecubicspline.hh>
#include <dumux/io/gnuplotinterface.hh>

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

    const auto test = [](auto f, auto df, const auto& testPoints, const auto& samplePoints, const std::string& prefix)
    {
        // create some test samples
        const auto ref = eval(f, testPoints);
        const auto refDeriv = eval(df, testPoints);

        // create the spline sample points
        const auto y = eval(f, samplePoints);

        // create the spline
        Dumux::MonotoneCubicSpline<double> spline(samplePoints, y);

        // evaluate spline and derivative
        const auto result = eval([&](const double x) { return spline.eval(x); }, testPoints);
        const auto resultDeriv = eval([&](const double x) { return spline.evalDerivative(x); }, testPoints);

        // compute largest difference
        auto diff = result; auto diffDeriv = result;
        std::transform(result.begin(), result.end(), ref.begin(), diff.begin(), [](auto a, auto b){ return std::abs(a-b); });
        std::transform(resultDeriv.begin(), resultDeriv.end(), refDeriv.begin(), diffDeriv.begin(), [](auto a, auto b){ return std::abs(a-b); });
        const auto maxNorm = std::accumulate(diff.begin(), diff.end(), diff[0], [](auto a, auto b){ return std::max(a, b); })
                             /std::abs(*std::max_element(ref.begin(), ref.end(), [](auto a, auto b){ return abs(a) < abs(b); }));
        const auto maxNormDeriv = std::accumulate(diffDeriv.begin(), diffDeriv.end(), diffDeriv[0], [](auto a, auto b){ return std::max(a, b); })
                                  /std::abs(*std::max_element(refDeriv.begin(), refDeriv.end(), [](auto a, auto b){ return abs(a) < abs(b); }));
        std::cout << "Maximum error: " << std::scientific << maxNorm << "\n";
        std::cout << "Maximum error in derivative: " << std::scientific << maxNormDeriv << "\n";

        if (maxNorm > 0.0008 || maxNormDeriv > 0.013)
            DUNE_THROW(Dune::Exception, "Maximum error in spline interpolation too large!");

        // test inverse by evaluating (x = f^-1(f(x))) for monotonically increasing function
        {
            const auto resultX = eval([&](double x){ return spline.evalInverse(spline.eval(x)); }, testPoints);
            auto diffInverse = resultX;
            std::transform(resultX.begin(), resultX.end(), testPoints.begin(), diffInverse.begin(), [](auto a, auto b){ return std::abs(a-b); });
            const auto maxNormInverse = std::accumulate(diffInverse.begin(), diffInverse.end(), diffInverse[0], [](auto a, auto b){ return std::max(a, b); })
                                 /std::abs(*std::max_element(testPoints.begin(), testPoints.end(), [](auto a, auto b){ return abs(a) < abs(b); }));


            std::cout << "Maximum error in identity using the inverse (mon. incr.): " << std::scientific << maxNormInverse << "\n";
            if (maxNormInverse > 1e-13)
                DUNE_THROW(Dune::Exception, "Maximum error in spline interpolation too large!");
        }
        // test inverse by evaluating (x = f^-1(f(x))) for monotonically decreasing function
        {
            auto reverseTest = testPoints;
            std::reverse(reverseTest.begin(), reverseTest.end());
            const auto resultX = eval([&](double x){ return spline.evalInverse(spline.eval(x)); }, reverseTest);
            auto diffInverse = resultX;
            std::transform(resultX.begin(), resultX.end(), reverseTest.begin(), diffInverse.begin(), [](auto a, auto b){ return std::abs(a-b); });
            const auto maxNormInverse = std::accumulate(diffInverse.begin(), diffInverse.end(), diffInverse[0], [](auto a, auto b){ return std::max(a, b); })
                                 /std::abs(*std::max_element(reverseTest.begin(), reverseTest.end(), [](auto a, auto b){ return abs(a) < abs(b); }));

            std::cout << "Maximum error in identity using the inverse (mon. decr.): " << std::scientific << maxNormInverse << "\n";
            if (maxNormInverse > 1e-13)
                DUNE_THROW(Dune::Exception, "Maximum error in spline interpolation too large!");
        }

        // plot with Gnuplot (plot a bit more so we can see the linear extension)
        const auto plotPoints = Dumux::linspace(-1.0, 5.0, 1000);
        const auto refPlot = eval(f, plotPoints);
        const auto refDerivPlot = eval(df, plotPoints);
        const auto resultPlot = eval([&](const double x) { return spline.eval(x); }, plotPoints);
        const auto resultDerivPlot = eval([&](const double x) { return spline.evalDerivative(x); }, plotPoints);

        Dumux::GnuplotInterface<double> gnuplot(/*persist=*/false);
        gnuplot.setOpenPlotWindow(false);
        gnuplot.addDataSetToPlot(plotPoints, refPlot, prefix + "exp_reference");
        gnuplot.addDataSetToPlot(plotPoints, refDerivPlot, prefix + "exp_reference_derivative");
        gnuplot.addDataSetToPlot(plotPoints, resultPlot, prefix + "monotspline");
        gnuplot.addDataSetToPlot(plotPoints, resultDerivPlot, prefix + "monotspline_derivative");
        gnuplot.plot(prefix + "monotspline");
    };

    {
        // we test the spline interpolation against a sample function
        // monotonically increasing function
        const auto f = [](double x){ return x*x*x; };
        const auto df = [](double x){ return 3*x*x; };

        const auto testPoints = Dumux::linspace(0.0, 4.0, 1000);
        const auto samplePoints = Dumux::linspace(0.0, 5.0, 10);
        test(f, df, testPoints, samplePoints, "x3in");

        const auto testPoints2 = Dumux::linspace(4.0, 0.0, 1000);
        const auto samplePoints2 = Dumux::linspace(5.0, 0.0, 10);
        test(f, df, testPoints2, samplePoints2, "x3de");
    }

    {
        // we test the spline interpolation against a sample function
        // monotonically decreasing function
        const auto f = [](double x){ return -x*x*x; };
        const auto df = [](double x){ return -3*x*x; };

        const auto testPoints = Dumux::linspace(0.0, 4.0, 1000);
        const auto samplePoints = Dumux::linspace(0.0, 5.0, 10);
        test(f, df, testPoints, samplePoints, "mx3in");

        const auto testPoints2 = Dumux::linspace(4.0, 0.0, 1000);
        const auto samplePoints2 = Dumux::linspace(5.0, 0.0, 10);
        test(f, df, testPoints2, samplePoints2, "mx3de");
    }

    return 0;
}
