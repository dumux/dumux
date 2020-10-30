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
 * \ingroup MaterialTests
 * \brief Test for two-phase material laws
 */

#include <config.h>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/splinemateriallaw.hh>

namespace Dumux {

template<class Function>
std::vector<double> eval(const Function& f, const std::vector<double>& x)
{
    auto y = x;
    std::transform(x.begin(), x.end(), y.begin(), [&](const double x) { return f(x); });
    return y;
}

template<class Function, class Orig, class Spline>
void runTest(const std::string& name, const Function& f,
             const Orig& orig, const Spline& spline,
             bool plot = true)
{
    std::cout << "-----------------------\n"
              << name << "\n"
              << "-----------------------\n";


    if (plot)
    {
        const auto [swMinPlot, swMaxPlot]
            = getParam<std::array<double, 2>>("Plot.Range", std::array<double, 2>{{0.1, 1.0}});

        const auto sw = linspace(swMinPlot, swMaxPlot, 1000);
        const auto value = eval([&](auto sw){ return f(sw, orig); }, sw);
        const auto valueSpline = eval([&](auto sw){ return f(sw, spline); }, sw);

        GnuplotInterface<double> gnuplot(false);
        gnuplot.addDataSetToPlot(sw, value, name + "-orig.dat", "title 'orig. curve'");
        gnuplot.addDataSetToPlot(sw, valueSpline, name + "-spline.dat", "title 'spline'");
        gnuplot.setOption("set xrange [0.0 : 1.0]");
        gnuplot.plot(name);
    }

    // speed test
    {
        const auto [swMinSpline, swMaxSpline]
            = getParam<std::array<double, 2>>("MaterialLaw.SplineSweInterval", std::array<double, 2>{{0.1, 1.0}});

        constexpr std::size_t testSamples = 1000000;
        const auto swTest = linspace(swMinSpline, swMaxSpline, testSamples);
        auto pcTest = swTest;

        Dune::Timer timer;
        auto resOrig = swTest;
        for (int i = 0; i < testSamples; ++i)
            resOrig[i] = f(swTest[i], orig);
        timer.stop();

        const auto vgTime = timer.elapsed();
        std::cout << "Unregularized law computed " << testSamples << " samples in " << vgTime << " seconds." << std::endl;

        timer.reset();
        timer.start();
        auto resSpline = swTest;
        for (int i = 0; i < testSamples; ++i)
            resSpline[i] = f(swTest[i], spline);
        timer.stop();

        const auto vgSplineTime = timer.elapsed();
        std::cout << "Spline law computed " << testSamples << " samples in " << vgSplineTime << " seconds." << std::endl;

        std::cout << "Speed-up factor ca. " << (vgTime/vgSplineTime) << "x (only in spline region)" << std::endl;

        auto error = resOrig;
        for (int i = 0; i < error.size(); ++i)
            error[i] = std::abs(error[i]-resSpline[i]);
        std::cout << "Maximum error: " << *std::max_element(error.begin(), error.end()) << std::endl;
    }
}

} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    Parameters::init(argc, argv);

    using MaterialLaw = FluidMatrix::VanGenuchtenNoReg<double>;
    MaterialLaw vg("MaterialLaw"); // read parameters from input file (group MaterialLaw)

    using MaterialLawSpline = FluidMatrix::SplineTwoPMaterialLaw<MaterialLaw>;
    MaterialLawSpline vgSpline("MaterialLaw"); // read parameters from input file (group MaterialLaw)

    const bool plot = getParam<bool>("Plot.EnablePlot");

    runTest("vg-pc", [](auto sw, const auto& law){ return law.pc(sw); }, vg, vgSpline, plot);
    runTest("vg-dpc", [](auto sw, const auto& law){ return law.dpc_dsw(sw); }, vg, vgSpline, plot);
    runTest("vg-krw", [](auto sw, const auto& law){ return law.krw(sw); }, vg, vgSpline, plot);
    runTest("vg-dkrw", [](auto sw, const auto& law){ return law.dkrw_dsw(sw); }, vg, vgSpline, plot);
    runTest("vg-krn", [](auto sw, const auto& law){ return law.krn(sw); }, vg, vgSpline, plot);
    runTest("vg-dkrn", [](auto sw, const auto& law){ return law.dkrn_dsw(sw); }, vg, vgSpline, plot);

    // inversions
    runTest("vg-sw-i", [](auto sw, const auto& law){ return law.sw(law.pc(sw)); }, vg, vgSpline, plot);
    runTest("vg-dsw-i", [](auto sw, const auto& law){ return law.dsw_dpc(law.pc(sw)); }, vg, vgSpline, plot);

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {

    std::cout << e << std::endl;
    return 1;
}
