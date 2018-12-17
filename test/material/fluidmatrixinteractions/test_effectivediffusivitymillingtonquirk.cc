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
 *
 * \brief Test for the Millington and Quirk effective diffusivity model
 */
#include <config.h>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/ploteffectivediffusivitymodel.hh>

#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    GnuplotInterface<double> gnuplot;
    gnuplot.setOpenPlotWindow(false);

    using EffDiffModel = DiffusivityMillingtonQuirk<double>;
    PlotEffectiveDiffusivityModel<double, EffDiffModel> plotEffectiveDiffusivityModel;
    const std::string fileName = "millingtonquirk_d_eff.dat";
    const double porosity = 0.3; // [-]
    plotEffectiveDiffusivityModel.adddeffcurve(gnuplot, porosity, 0.0, 1.0, fileName);

    gnuplot.plot("d_eff");

    return 0;
}
