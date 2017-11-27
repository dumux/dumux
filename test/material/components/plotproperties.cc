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
 * \brief Plot properties of components and fluids
 */

#include "config.h"
#include <array>
#include <cstring>
#include <limits>
#include <vector>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/benzene.hh>
#include <dumux/material/components/brine.hh>
#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/material/components/h2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/heavyoil.hh>
#include <dumux/material/components/lnapl.hh>
#include <dumux/material/components/mesitylene.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/simpleco2.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/xylene.hh>

using namespace std;
////////////////////////
// the main function
////////////////////////
int main(int argc, char *argv[])
{
    bool openPlotWindow = true;
    if (argc == 2 && (strcmp(argv[1], "0") || strcmp(argv[1], "false") || strcmp(argv[1], "False")))
        openPlotWindow = false;

    double pressure = 1e5;
    double TMin = 273.15;
    double TMax = 323.15;
    double TRange = TMax - TMin;
    const unsigned int numIntervals = 100;
    vector<double> T(numIntervals + 1);
    for (int i = 0; i <= numIntervals; i++)
    {
        T[i] = TMin + TRange * double(i) /double(numIntervals);
    }

    // components
    const unsigned int liquidPhaseIdx = 0;
    const unsigned int gasPhaseIdx = 1;
    using Component = COMPONENT;

    const unsigned int numPhases = 2;
    array<string, numPhases> phaseNames;
    phaseNames[liquidPhaseIdx] = "liquid";
    phaseNames[gasPhaseIdx] = "gas";

    const unsigned int numProperties = 5;
    array<string, numProperties> propertyNames;
    array<string, numProperties> propertyUnits;
    unsigned int densityIdx = 0;
    propertyNames[densityIdx] = "density";
    propertyUnits[densityIdx] = "[kg/m^3]";
    unsigned int enthalpyIdx = 1;
    propertyNames[enthalpyIdx] = "enthalpy";
    propertyUnits[enthalpyIdx] = "[J/(kg)]";
    unsigned int heatCapacityIdx = 2;
    propertyNames[heatCapacityIdx] = "heatCapacity";
    propertyUnits[heatCapacityIdx] = "[J/(kg*K)]";
    unsigned int viscosityIdx = 3;
    propertyNames[viscosityIdx] = "viscosity";
    propertyUnits[viscosityIdx] = "[Pa*s]";
    unsigned int thermalConductivityIdx = 4;
    propertyNames[thermalConductivityIdx] = "thermalConductivity";
    propertyUnits[thermalConductivityIdx] = "[W/(m*K)]";
    array<array<vector<double>, numProperties>, numPhases> property;
    array<array<bool, numProperties>, numPhases> propertyAvailable;
    array<array<array<double, 2>, numProperties>, numPhases> propertyMinMax;

    // get values from component functions
    for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        for (unsigned int propertyIdx = 0; propertyIdx < numProperties; ++propertyIdx)
        {
            propertyAvailable[phaseIdx][propertyIdx] = true;
            propertyMinMax[phaseIdx][propertyIdx][0] = std::numeric_limits<double>::max();
            propertyMinMax[phaseIdx][propertyIdx][1] = std::numeric_limits<double>::min();
            property[phaseIdx][propertyIdx].resize(numIntervals+1);
        }

        for (int i = 0; i <= numIntervals; i++)
        {
            if (phaseIdx == liquidPhaseIdx)
            {
                try { property[phaseIdx][densityIdx][i] = Component::liquidDensity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][densityIdx] = false; }
                try { property[phaseIdx][enthalpyIdx][i] = Component::liquidEnthalpy(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][enthalpyIdx] = false; }
                try { property[phaseIdx][heatCapacityIdx][i] = Component::liquidHeatCapacity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][heatCapacityIdx] = false; }
                try { property[phaseIdx][viscosityIdx][i] = Component::liquidViscosity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][viscosityIdx] = false; }
                try { property[phaseIdx][thermalConductivityIdx][i] = Component::liquidThermalConductivity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][thermalConductivityIdx] = false; }
            }
            if (phaseIdx == gasPhaseIdx)
            {
                try { property[phaseIdx][densityIdx][i] = Component::gasDensity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][densityIdx] = false; }
                try { property[phaseIdx][enthalpyIdx][i] = Component::gasEnthalpy(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][enthalpyIdx] = false; }
                try { property[phaseIdx][heatCapacityIdx][i] = Component::gasHeatCapacity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][heatCapacityIdx] = false; }
                try { property[phaseIdx][viscosityIdx][i] = Component::gasViscosity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][viscosityIdx] = false; }
                try { property[phaseIdx][thermalConductivityIdx][i] = Component::gasThermalConductivity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][thermalConductivityIdx] = false; }
            }

            for (unsigned int propertyIdx = 0; propertyIdx < numProperties; ++propertyIdx)
            {
                propertyMinMax[phaseIdx][propertyIdx][0] = std::min(propertyMinMax[phaseIdx][propertyIdx][0], property[phaseIdx][propertyIdx][i]);
                propertyMinMax[phaseIdx][propertyIdx][1] = std::max(propertyMinMax[phaseIdx][propertyIdx][1], property[phaseIdx][propertyIdx][i]);
            }
        }
    }

    // plot densities
    for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        for (unsigned int propertyIdx = 0; propertyIdx < numProperties; ++propertyIdx)
        {
            if (!propertyAvailable[phaseIdx][propertyIdx])
                continue;

            Dumux::GnuplotInterface<double> gnuplot(true);
            gnuplot.setOpenPlotWindow(openPlotWindow);
            gnuplot.setCreateImage(true);
            gnuplot.setXRange(TMin, TMax);
            gnuplot.setYRange(propertyMinMax[phaseIdx][propertyIdx][0]*0.999, propertyMinMax[phaseIdx][propertyIdx][1]*1.001);
            gnuplot.setXlabel("temperature [K]");
            gnuplot.setYlabel(phaseNames[phaseIdx] + " " + propertyNames[propertyIdx] + " " + propertyUnits[propertyIdx]);
            gnuplot.setDatafileSeparator(',');
            gnuplot.addDataSetToPlot(T, property[phaseIdx][propertyIdx], Component::name() + "_" + phaseNames[phaseIdx] + "_" + propertyNames[propertyIdx] + ".csv");
            gnuplot.plot(Component::name() + "_" + phaseNames[phaseIdx] + "_" + propertyNames[propertyIdx]);
        }
    }
}
