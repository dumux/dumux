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

#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/datasplinemateriallaw.hh>

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
}

} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    Parameters::init(argc, argv);

    using MaterialLaw = FluidMatrix::LinearMaterialDefault<double>;
    MaterialLaw lin("MaterialLaw"); // read parameters from input file (group MaterialLaw)

    using MaterialLawSpline = FluidMatrix::DataSplineTwoPMaterialLaw<double>;
    MaterialLawSpline linSpline("MaterialLaw"); // read parameters from input file (group MaterialLaw)

    const bool plot = getParam<bool>("Plot.EnablePlot");

    runTest("lin-pc", [](auto sw, const auto& law){ return law.pc(sw); }, lin, linSpline, plot);
    runTest("lin-dpc", [](auto sw, const auto& law){ return law.dpc_dsw(sw); }, lin, linSpline, plot);
    runTest("lin-krw", [](auto sw, const auto& law){ return law.krw(sw); }, lin, linSpline, plot);
    runTest("lin-dkrw", [](auto sw, const auto& law){ return law.dkrw_dsw(sw); }, lin, linSpline, plot);
    runTest("lin-krn", [](auto sw, const auto& law){ return law.krn(sw); }, lin, linSpline, plot);
    runTest("lin-dkrn", [](auto sw, const auto& law){ return law.dkrn_dsw(sw); }, lin, linSpline, plot);

    // inversions
    runTest("lin-sw-i", [](auto sw, const auto& law){ return law.sw(law.pc(sw)); }, lin, linSpline, plot);
    runTest("lin-dsw-i", [](auto sw, const auto& law){ return law.dsw_dpc(law.pc(sw)); }, lin, linSpline, plot);

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {

    std::cout << e << std::endl;
    return 1;
}
