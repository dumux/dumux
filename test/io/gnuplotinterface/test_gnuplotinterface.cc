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
 * \brief Simple test for the gnuplot interface, which plots a parabola
 */
#include "config.h"
#include <vector>
#include <dumux/io/gnuplotinterface.hh>

////////////////////////
// the main function
////////////////////////
int main()
{
    Dumux::GnuplotInterface<double> gnuplot;

    unsigned int numIntervals = 100;
    std::vector<double> x(101);
    std::vector<double> y(101);
    double xMin = -5.0;
    double xMax = 5.0;
    double xRange = xMax - xMin;
    double yMin = 1e100;
    double yMax = -1e100;

    for (int i = 0; i <= numIntervals; i++)
    {
        x[i] = xMin + xRange * double(i) /double(numIntervals);
        y[i] = x[i]*x[i];
        yMin = std::min(yMin, y[i]);
        yMax = std::max(yMax, y[i]);
    }

    unsigned int windowNumber = 0;
    gnuplot.reset(windowNumber);
    gnuplot.setXRange(0, 5, windowNumber);
    gnuplot.setYRange(yMin, yMax, windowNumber);
    gnuplot.setXlabel("x [-]", windowNumber);
    gnuplot.setYlabel("f(x) [-]", windowNumber);
    gnuplot.addDataSetToPlot(x, y, "f(x)=x^2", windowNumber);
    gnuplot.addFunctionToPlot("x**3", "f(x)=x^3", windowNumber);
    gnuplot.addFileToPlot("0_f(x)=x^2.dat", "0_f(x)=x^2.dat", windowNumber);
    gnuplot.plot("plot", windowNumber, true);
    exit(0);
}
