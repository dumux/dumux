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
#if AIR
    typedef Dumux::Air<double> Component;
#elif BENZENE
    typedef Dumux::Benzene<double> Component;
#elif BRINE
    typedef Dumux::Brine<double> Component;
#elif METHANE
    typedef Dumux::CH4<double> Component;
#elif CARBONDIOXIDE
    typedef Dumux::CO2<double> Component;
#elif DNAPL_TCE
    typedef Dumux::DNAPL<double> Component;
#elif HYDROGEN
    typedef Dumux::H2<double> Component;
#elif WATER
    typedef Dumux::H2O<double> Component;
#elif HEAVYOIL
    typedef Dumux::HeavyOil<double> Component;
#elif LNAPL_OIL
    typedef Dumux::LNAPL<double> Component;
#elif MESITYLENE
    typedef Dumux::Mesitylene<double> Component;
#elif NITROGEN
    typedef Dumux::N2<double> Component;
#elif OXYGEN
    typedef Dumux::O2<double> Component;
#elif SIMPLE_CARBONDIOXIDE
    typedef Dumux::SimpleCO2<double> Component;
#elif SIMPLE_WATER
    typedef Dumux::SimpleH2O<double> Component;
#elif XYLENE
    typedef Dumux::Xylene<double> Component;
#endif

    const unsigned int numPhases = 2;
    array<string, numPhases> phaseNames;
    phaseNames[liquidPhaseIdx] = "liquid";
    phaseNames[gasPhaseIdx] = "gas";

    const unsigned int numProperties = 4;
    array<string, numProperties> propertyNames;
    array<string, numProperties> propertyUnits;
    propertyNames[0] = "density";
    propertyUnits[0] = "[kg/m^3]";
    propertyNames[1] = "heatCapacity";
    propertyUnits[1] = "[J/(kg*K)]";
    propertyNames[2] = "viscosity";
    propertyUnits[2] = "[Pa*s]";
    propertyNames[3] = "thermalConductivity";
    propertyUnits[3] = "[W/(m*K)]";
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
                try { property[phaseIdx][0][i] = Component::liquidDensity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][0] = false; }
                try { property[phaseIdx][1][i] = Component::liquidHeatCapacity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][1] = false; }
                try { property[phaseIdx][2][i] = Component::liquidViscosity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][2] = false; }
                try { property[phaseIdx][3][i] = Component::liquidThermalConductivity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][3] = false; }
            }
            if (phaseIdx == gasPhaseIdx)
            {
                try { property[phaseIdx][0][i] = Component::gasDensity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][0] = false; }
                try { property[phaseIdx][1][i] = Component::gasHeatCapacity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][1] = false; }
                try { property[phaseIdx][2][i] = Component::gasViscosity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][2] = false; }
                try { property[phaseIdx][3][i] = Component::gasThermalConductivity(T[i], pressure); }
                catch (Dune::NotImplemented &e) { propertyAvailable[phaseIdx][3] = false; }
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
