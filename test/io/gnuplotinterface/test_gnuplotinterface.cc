// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Simple test for the gnuplot interface, which plots a parabola
 */
#include <config.h>
#include <vector>
#include <dumux/io/gnuplotinterface.hh>

////////////////////////
// the main function
////////////////////////
int main()
{
    Dumux::GnuplotInterface<double> gnuplot(false);
    gnuplot.setOpenPlotWindow(false);

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
        using std::max;
        using std::min;
        yMin = min(yMin, y[i]);
        yMax = max(yMax, y[i]);
    }

    gnuplot.setOutputDirectory("output");
    gnuplot.setXRange(0, 5);
    gnuplot.setYRange(yMin, yMax);
    gnuplot.setXlabel("x [-]");
    gnuplot.setYlabel("f(x) [-]");
    gnuplot.setDatafileSeparator(',');
    gnuplot.setOption("set arrow from 0,0 to 2,20 head filled lc rgb 'dark-gray'");
    gnuplot.setOption("set label 'arrow' at 1,15 center tc rgb 'dark-gray'");
    gnuplot.addDataSetToPlot(x, y, "dataSet.csv", "every 5 w lp ps 2");
    gnuplot.addFunctionToPlot("x**3", "title 'function_f(x)=x^3'");
    gnuplot.addFileToPlot("dataSet.csv");
    gnuplot.plot("plot");
    exit(0);
}
