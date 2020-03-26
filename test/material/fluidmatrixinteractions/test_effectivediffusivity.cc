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
 * \brief Test for the Millington and Quirk effective diffusivity model.
 */

#include <config.h>

#include <dumux/common/parameters.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/ploteffectivediffusivitymodel.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    Parameters::init(argc, argv);

    GnuplotInterface<double> gnuplot;
    gnuplot.setOpenPlotWindow(false);
    gnuplot.setXlabel("phase saturation [-]");
    gnuplot.setYlabel("effective diffusion/molecular diffusion [-]");

    const double porosity = 0.3; // [-]
    const double diffCoeff = 1.0;

    using ConstantEffDiffModel = DiffusivityConstantTortuosity<double>;
    PlotEffectiveDiffusivityModel<double, ConstantEffDiffModel> plotConstantTortuosityDiffusivityModel;
    const std::string fileNameConstant = "constant_d_eff.dat";
    plotConstantTortuosityDiffusivityModel.adddeffcurve(gnuplot, porosity, diffCoeff, 0.0, 1.0, fileNameConstant);

    using MillingtonQuirkEffDiffModel = DiffusivityMillingtonQuirk<double>;
    PlotEffectiveDiffusivityModel<double, MillingtonQuirkEffDiffModel> plotMillingtonQuirkDiffusivityModel;
    const std::string fileNameMQ = "millingtonquirk_d_eff.dat";
    plotMillingtonQuirkDiffusivityModel.adddeffcurve(gnuplot, porosity, diffCoeff, 0.0, 1.0, fileNameMQ);

    gnuplot.plot("d_eff");

    // Test against one value
    double millingtonQuirk_d_eff = plotMillingtonQuirkDiffusivityModel.getEffectiveDiffusionCoefficient(1.0,1.0,1.0);
    double constantTortuoisty_d_eff = plotConstantTortuosityDiffusivityModel.getEffectiveDiffusionCoefficient(1.0,1.0,1.0);
    if (constantTortuoisty_d_eff != 0.3)
        DUNE_THROW(Dune::Exception, "The constant tortuosity method for calculating the effective diffusion coefficient is not correct");
    if (millingtonQuirk_d_eff != 1.0)
        DUNE_THROW(Dune::Exception, "The Millington Quirk method for calculating the effective diffusion coefficient is not correct");

    return 0;
}
