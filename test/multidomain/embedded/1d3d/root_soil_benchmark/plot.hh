// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \author Timo Koch <timokoch@uio.no>
 * \brief Plotting tools for root soil benchmark
 */
#ifndef DUMUX_TEST_ROOT_SOIL_BENCHMARK_PLOT_HH
#define DUMUX_TEST_ROOT_SOIL_BENCHMARK_PLOT_HH

#include <vector>
#include <algorithm>
#include <tuple>
#include <utility>
#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/io/gnuplotinterface.hh>

namespace Dumux::RootSoil {

/*!
 * \brief plot the transpiration curves
 * \tparam Problem the problem type
 */
template<class Problem>
class TranspirationPlot
{
public:
    TranspirationPlot(std::shared_ptr<const Problem> problem)
    : problem_(problem)
    {
        const std::string outputDir = getParam<std::string>("Output.GnuplotOutputDirectory", "");
        filenamePrefix_ = getParam<std::string>("Output.GnuplotOutputFilenamePrefix", "");
        const bool openGnuPlot = getParam<bool>("Output.GnuplotShow", true);
        gnuplot_.setOpenPlotWindow(openGnuPlot);
        if (outputDir != "")
            gnuplot_.setOutputDirectory(outputDir);
    }

    /*!
     * \brief add a data point to the transpiration plot
     * \param curSol the current solution vector (root)
     * \param curGridVars the current grid variables (root)
     * \param t the time at the end of the current time step
     * \param dt the time step size of the current time step
     */
    template<class SolutionVector, class GridVariables>
    void addDataPoint(const SolutionVector& curSol, const GridVariables& curGridVars, double t, double dt)
    {
        gnuplot_.resetPlot();

        gnuplot_.setXlabel("time [day]");
        gnuplot_.setYlabel("transpiration rate [g/day]");
        gnuplot_.setOption("set y2label \"cumulative transpiration [g]\"");

        const auto conversion = 86400.0*1000.0; // convert to g/day
        const auto tAct = problem_->computeActualTranspirationRate(curSol, curGridVars)*conversion;

        tActPlot_.push_back(tAct);
        if (tCumulPlot_.empty())
            tCumulPlot_.push_back(tAct*dt/86400.0);
        else
            tCumulPlot_.push_back(tCumulPlot_.back() + tAct*dt/86400.0);
        timePlot_.push_back(t/86400.0);
        gnuplot_.addDataSetToPlot(timePlot_, tActPlot_, filenamePrefix_ + "actualtranspiration.out", "with lines axes x1y1 lw 3");
        gnuplot_.addDataSetToPlot(timePlot_, tCumulPlot_, filenamePrefix_ + "cumulativetranspiration.out", "with lines axes x1y2 lw 3");

        if (problem_->bcType() != Problem::BCType::constCollarPressure)
        {
            const auto tPot = problem_->potentialTranspirationRate()*conversion;
            tPotPlot_.push_back(tPot);
            gnuplot_.addDataSetToPlot(timePlot_, tPotPlot_, filenamePrefix_ + "potentialtranspiration.out", "with lines axes x1y1 lw 2 lc rgb 'black'");
        }

        gnuplot_.setOption("set ytics nomirror");
        gnuplot_.setOption("set y2tics");

        gnuplot_.setOption("set autoscale x");
        gnuplot_.setOption("set autoscale y");
        gnuplot_.setOption("set autoscale y2");

        gnuplot_.setOption("set title \"Plant transpiration\"");
        gnuplot_.plot(filenamePrefix_ + "transpiration");
    }

private:
    std::shared_ptr<const Problem> problem_;
    GnuplotInterface<double> gnuplot_;
    std::string filenamePrefix_;

    // the cached data vectors (they are updated on each call to addDataPoint)
    std::vector<double> tPotPlot_, tActPlot_, tCumulPlot_, timePlot_;
};

} // end namespace Dumux::RootSoil

#endif
