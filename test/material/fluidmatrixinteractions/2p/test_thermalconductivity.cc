// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief Test for the Johansen thermal conductivity law.
 */

#include <config.h>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/plotthermalconductivitymodel.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/johansen.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>

#include <dumux/material/fluidsystems/h2on2.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    using Scalar = double;
    using FluidSystem = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;

    GnuplotInterface<double> gnuplot;
    gnuplot.setOpenPlotWindow(false);
    gnuplot.setXlabel("wetting phase saturation [-]");
    gnuplot.setYlabel("effective thermal conductivity [W/(m K)]");
    gnuplot.setOption("set key right bottom reverse samplen 1");

    const Scalar porosity = 0.3; // [-]
    const Scalar rhoSolid = 2700.0; // kg/m^3
    const Scalar lambdaSolid = 2.8; // W/(m K)

    using JohansenThermCondModel = ThermalConductivityJohansen<Scalar>;
    PlotThermalConductivityModel<Scalar, JohansenThermCondModel, FluidSystem> plotJohansenThermalConductivityModel(293.15, 1e5);
    const std::string fileNameJohansen = "johansen_lambda_eff.dat";
    plotJohansenThermalConductivityModel.addlambdaeffcurve(gnuplot, porosity, rhoSolid, lambdaSolid, 0.0, 1.0, fileNameJohansen);

    using SomertonThermCondModel = ThermalConductivitySomerton<Scalar>;
    PlotThermalConductivityModel<Scalar, SomertonThermCondModel, FluidSystem> plotSomertonThermalConductivityModel(293.15, 1e5);
    const std::string fileNameSomerton = "somerton_lambda_eff.dat";
    plotSomertonThermalConductivityModel.addlambdaeffcurve(gnuplot, porosity, rhoSolid, lambdaSolid, 0.0, 1.0, fileNameSomerton);

    gnuplot.plot("lambda_eff");

    return 0;
}
