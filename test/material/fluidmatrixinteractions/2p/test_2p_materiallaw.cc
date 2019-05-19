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
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/common/parameters.hh>

#include <dumux/material/fluidmatrixinteractions/2pnew/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2pnew/materiallaw.hh>

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

int main(int argc, char** argv) try
{
    using namespace Dumux;

    Parameters::init(argc, argv);

    using MaterialLaw = FluidMatrix::TwoPMaterialLaw<double, FluidMatrix::VanGenuchten, FluidMatrix::NoTwoPRegularization<double>>;
    auto vg = std::make_shared<MaterialLaw>("MaterialLaw"); // read parameters from input file (group MaterialLaw)

    const auto swMinPlot = getParam<double>("Plot.SwMin", 0.1);
    const auto swMaxPlot = getParam<double>("Plot.SwMax", 0.9);
    const auto sw = linspace(swMinPlot,swMaxPlot, 1000);
    const auto pc = eval([&](auto sw){ return vg->pc(sw); }, sw);

    GnuplotInterface<double> gnuplot(true);
    gnuplot.addDataSetToPlot(sw, pc, "pcsw.dat", "title 'van Genuchten pc-sw curve'");
    gnuplot.setOption("set xrange [0.0 : 1.0]");
    gnuplot.plot("vangenuchten");

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {

    std::cout << e << std::endl;
    return 1;
}
