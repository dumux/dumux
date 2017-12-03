// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Test for the Somerton thermal conductivity law
 */
 #include <config.h>

 #include <dumux/io/gnuplotinterface.hh>
 #include <dumux/io/plotthermalconductivitymodel.hh>

 #include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
 #include <dumux/material/fluidsystems/h2on2.hh>

 namespace Dumux {
 namespace Properties {
 NEW_TYPE_TAG(TestTypeTag);
 SET_TYPE_PROP(TestTypeTag, Scalar, double);
 SET_TYPE_PROP(TestTypeTag, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);
 SET_TYPE_PROP(TestTypeTag, Indices, SomertonIndices);
 SET_TYPE_PROP(TestTypeTag, ThermalConductivityModel, ThermalConductivitySomerton<typename GET_PROP_TYPE(TypeTag, Scalar)>);
 } // end namespace Properties
 } // end namespace Dumux

 int main(int argc, char** argv)
 {
     using namespace Dumux;
     using TypeTag = TTAG(TestTypeTag);

     GnuplotInterface<double> gnuplot;
     gnuplot.setOpenPlotWindow(false);

     PlotThermalConductivityModel<TypeTag> plotThermalConductivityModel(293.15, 1e5);
     const std::string fileName = "somerton_lambda_eff.dat";
     const double porosity = 0.3; // [-]
     const double rhoSolid = 2700.0; // kg/m^3
     const double lambdaSolid = 2.8; // W/(m K)
     plotThermalConductivityModel.addlambdaeffcurve(gnuplot, porosity, rhoSolid, lambdaSolid, 0.0, 1.0, fileName);

     gnuplot.plot("lambda_eff");

     return 0;
 }
