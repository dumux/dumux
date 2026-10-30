// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_HELPER_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_HELPER_HH

/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief A helper class that determines the volume flux leaving the pore network in the flow
 *        direction, converts it into an apparent permeability, and collects/plots the
 *        resulting relative-permeability curves versus saturation.
 *
 * \note The pc-Sw curve itself (common to both phases) is plotted once per cycle directly from
 *       the quasi-static two-phase sweep, see PNMProblemStatic::plot() in problem_static.hh --
 *       it is not duplicated here.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <functional>
#include <array>
#include <vector>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup PoreNetworkModels
 * \brief Collects, writes and plots apparent permeability vs. saturation data gathered over
 *        the course of an upscaling run.
 *
 * One instance is shared across both single-phase (wetting/non-wetting) upscaling solves, so
 * that both end up in the same plot; see main.cc.
 */
template<class Scalar>
class UpscalingHelper
{
public:
    //! Records the total pore volume of the network (needed for the porosity report).
    void setTotalPoreVolume(Scalar totalPoreVolume)
    { totalPoreVolume_ = totalPoreVolume; }

    /*!
     * \brief Records the row index at which the drainage leg ends and the imbibition leg
     *        begins in the collected data (i.e. the number of pc steps of the drainage leg),
     *        so the two legs can later be plotted as separate, distinctly styled curves.
     */
    void setApexStep(std::size_t apexStep)
    { apexStep_ = apexStep; }

    /*!
     * \brief Clears the collected data of the given direction.
     *
     * Since one UpscalingHelper instance is shared across the wetting- and non-wetting-phase
     * upscaling runs (see main.cc) so that both end up in the same plot/porosity report, this
     * must be called before starting a new phase's step loop for a given direction -- otherwise
     * that phase's data would be appended to, rather than replace, the previous phase's.
     */
    void resetDataPoints(std::size_t dirIdx)
    {
        effectivePermeability_[dirIdx].clear();
        averageSaturation_[dirIdx].clear();
    }

    /*!
     * \brief Evaluates the mass flux leaving the network in the direction of the applied
     *        pressure gradient, converts it into an apparent permeability, and stores it
     *        together with the corresponding saturation.
     */
    template<class Problem, class StaticProperties>
    void setDataPoints(const Problem& problem, const Scalar totalMassFlux, const StaticProperties& staticProperties)
    {
        // only the permeability in flow direction (set via Problem.Directions, x by default)
        // is upscaled; side lengths of the other directions are needed for the outflow area only
        auto sideLengths = problem.sideLengths();
        const auto direction = problem.direction();

        const auto density = problem.density();
        const auto dynamicViscosity = problem.dynamicViscosity();
        const auto pressureGradient = problem.pressureGradient();

        // convert mass to volume flux
        const auto volumeFlux = totalMassFlux / density;

        // outflow area of a cylindrical representative elementary volume perpendicular to
        // the flow direction
        sideLengths[direction] = 1.0;
        const auto outflowArea = M_PI/4 * std::accumulate(sideLengths.begin(), sideLengths.end(), Scalar(1.0), std::multiplies<Scalar>());
        const auto apparentVelocity = volumeFlux / outflowArea;
        const auto apparentPermeability = apparentVelocity / pressureGradient * dynamicViscosity;

        effectivePermeability_[direction].push_back(apparentPermeability);
        averageSaturation_[direction].push_back(staticProperties.averageSaturation);
    }

    //! Computes and stores the porosity of the network (pore volume over REV volume).
    template<class Problem>
    void networkPorosity(const Problem& problem)
    {
        const auto& sideLengths = problem.sideLengths();
        const auto networkTotalVolume = M_PI/4 * std::accumulate(sideLengths.begin(), sideLengths.end(), Scalar(1.0), std::multiplies<Scalar>());
        networkPorosity_ = totalPoreVolume_ / networkTotalVolume;
    }

    //! Reports the previously computed network porosity to the terminal.
    void reportPorosity() const
    { std::cout << "the porosity of the network is: " << networkPorosity_ << std::endl; }

    // We determine the domain side lengths by using the bounding box of the network.
    template<class GridGeometry>
    auto getSideLengths(const GridGeometry& gridGeometry) const
    {
        using GlobalPosition = typename GridGeometry::GlobalCoordinate;
        GlobalPosition result(0.0);

        std::cout << "Automatically determining side lengths of REV based on bounding box of pore network" << std::endl;
        for (int dimIdx = 0; dimIdx < GridGeometry::GridView::dimensionworld; ++dimIdx)
            result[dimIdx] = gridGeometry.bBoxMax()[dimIdx] - gridGeometry.bBoxMin()[dimIdx];

        return result;
    }

    //! Plots the collected effective permeability vs. saturation data.
    void plot() const
    { plotEffPermeabilityVsSw_(); }

    /*!
     * \brief Writes the collected effective-permeability data of the given direction and phase
     *        to disk.
     * \param dirIdx The spatial direction (0:x, 1:y, 2:z)
     * \param phaseIdx The phase (0: wetting, 1: non-wetting)
     * \param namePrefix A prefix (e.g. the run name) used to keep output files of different
     *                   runs from overwriting each other
     */
    void writePlotDataToFile(std::size_t dirIdx, int phaseIdx, const std::string& namePrefix)
    {
        const std::string fileName = namePrefix + "_" + dirName_[dirIdx] + "-dir-" + phaseName_[phaseIdx] + "-EffPermeabilityVsAverageSaturation.dat";
        std::ofstream logfile(fileName);
        for (std::size_t i = 0; i < effectivePermeability_[dirIdx].size(); ++i)
            logfile << averageSaturation_[dirIdx][i] << " " << effectivePermeability_[dirIdx][i] << std::endl;

        filesToPlot_.push_back(fileName);
    }

private:
    /*!
     * \brief Adds the drainage and imbibition legs of one data file to a gnuplot plot as two
     *        separate, distinctly styled and labeled curves (solid circles for drainage,
     *        dashed squares for imbibition; color picked per curve index so different files,
     *        e.g. different phases, remain visually distinguishable).
     */
    template<class Gnuplot>
    void addHysteresisLegsToPlot_(Gnuplot& gnuplot, const std::string& fileName, const std::string& curveLabel, int colorIndex) const
    {
        const std::string color = colorIndex % 2 == 0 ? "#1f77b4" : "#d62728";
        gnuplot.addFileToPlot(fileName, "every ::0::" + std::to_string(apexStep_) +
                              " u 1:2 title \"" + curveLabel + " (drainage)\" with linespoints lc rgb '" + color + "' pt 7 ps 0.6");
        gnuplot.addFileToPlot(fileName, "every ::" + std::to_string(apexStep_) +
                              " u 1:2 title \"" + curveLabel + " (imbibition)\" with linespoints lc rgb '" + color + "' dt 2 pt 5 ps 0.6");
    }

    void plotEffPermeabilityVsSw_() const
    {
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        gnuplot.resetAll();

        int colorIndex = 0;
        for (const auto& fileName : filesToPlot_)
        {
            const bool isWetting = fileName.find(phaseName_[0]) != std::string::npos;
            addHysteresisLegsToPlot_(gnuplot, fileName, isWetting ? "Wetting" : "Non-wetting", colorIndex++);
        }

        gnuplot.setOption("set title \"Effective permeability vs. Sw\"\n");
        gnuplot.setXlabel("Sw [-]");
        gnuplot.setYlabel("Effective Permeability [m2]");
        gnuplot.plot("Effective Permeability vs. Sw");
    }

    std::array<std::vector<Scalar>, 3> effectivePermeability_;
    std::array<std::vector<Scalar>, 3> averageSaturation_;
    const std::array<std::string, 3> dirName_ = {"X", "Y", "Z"};
    const std::array<std::string, 2> phaseName_{"Wetting", "Non-wetting"};
    std::vector<std::string> filesToPlot_;
    std::size_t apexStep_ = 0;
    Scalar totalPoreVolume_;
    Scalar networkPorosity_;
};

} // end namespace Dumux
#endif
