// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_HELPER_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_HELPER_HH

// ## Upscaling helper struct (`upscalinghelper.hh`)
//
// This file contains the __upscaling helper struct__ which considers the volume flux leaving
// the pore network in flow direction in order to find the upscaled Darcy permeability.
// [[content]]
#include <iostream>
#include <ostream>
#include <iomanip>
#include <numeric>
#include <functional>
#include <filesystem>

#include <dune/common/float_cmp.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/format.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

template<class Scalar>
class UpscalingHelper
{
public:
    // ### Set data points to calculate intrinsic permeability and Forchheimer coefficient
    // This function first evaluates the mass flux leaving the network in the direction of the applied pressure gradient.
    // Afterwards, the mass flux is converted into an volume flux which is used to calculate the apparent velocity.
    // Then apparent permeability of the network is computed and stored for furthure calculations.
    // [[codeblock]]
    template <class Problem, class StaticProperties>
    void setDataPoints(const Problem &problem, const Scalar totalMassFlux, const StaticProperties staticProperties)
    {
        // get the domain side lengths from the problem
        auto sideLengths = problem.sideLengths();


        // get the applied pressure gradient
        const auto pressureGradient = problem.pressureGradient();

        // get the fluid properties
        const auto density = problem.density();
        const auto dynamicViscosity = problem.dynamicViscosity();

        // convert mass to volume flux
        const auto volumeFlux = totalMassFlux / density;
        ;

        // calculate apparent velocity
        sideLengths[problem.direction()] = 1.0;
        const auto outflowArea = std::accumulate(sideLengths.begin(), sideLengths.end(), 1.0, std::multiplies<Scalar>());
        const auto vApparent = volumeFlux / outflowArea;

        // compute apparent permeability
        const auto KApparent = vApparent / pressureGradient * dynamicViscosity;

        // store the required data for further calculations
        EffectivePermeability_[problem.direction()].push_back(KApparent);
        averageSaturation_[problem.direction()].push_back(staticProperties.averageSaturation);
        pc_[problem.direction()].push_back(staticProperties.pcGlobal);
    }
    // [[/codeblock]]

    // ### Determine the domain's side lengths
    //
    // We determine the domain side length by using the bounding box of the network
    // [[codeblock]]
    template<class GridGeometry>
    auto getSideLengths(const GridGeometry& gridGeometry)
    {
        using GlobalPosition = typename GridGeometry::GlobalCoordinate;
        GlobalPosition result(0.0);

        std::cout << "Automatically determining side lengths of REV based on bounding box of pore network" << std::endl;
        for (int dimIdx = 0; dimIdx < GridGeometry::GridView::dimensionworld; ++dimIdx)
            result[dimIdx] = gridGeometry.bBoxMax()[dimIdx] - gridGeometry.bBoxMin()[dimIdx];

        return result;
    }
    // [[/codeblock]]

    // ### Plot the data using Gnuplot
    //
    // [[codeblock]]
    static void plot()
    {
        // // plot permeability ratio vs. Forchheimer number
        plotRelPermeabilityVsSw_();
        plotCapillaryPressureVsSw_();
    }
    // [[/codeblock]]
    //
    // ### Save the relevant data for plot
    // [[codeblock]]
    void writePlotDataToFile(std::size_t dirIdx, int phaseIdx)
    {
        // // write permeability ratio vs. Forchheimer number
        // writePermeabilityratioVsForchheimerNumber_(dirIdx);

        // // write inverse of apparent permability vs. rho v / mu
        // writeInversePrmeabilityVsInertiaToViscousRatio_(dirIdx);

        writeEffPermeabilityVsAverageSaturation_(dirIdx, phaseIdx );
        writeCapillaryPressureVsAverageSaturation_(dirIdx);
    }
    // [[/codeblock]]
    //
    // ### Report the upscaled data
    // [[codeblock]]
    void report(bool isCreepingFlow)
    {
        // Report the results for each direction
        // for (const auto dirIdx : directions_)
        // {
        //     std::cout << Fmt::format("\n{:#>{}}\n\n", "", 40)
        //               << Fmt::format("{}-direction:\n", dirName_[dirIdx])
        //               << Fmt::format("-- Darcy (intrinsic) permeability = {:.3e} m^2\n", EffectivePermeability_[dirIdx]);

        //     std::cout << Fmt::format("\n{:#>{}}\n", "", 40) << std::endl;
        // }
    }
    // [[/codeblock]]
    //
    // ### Compare with reference data provided in input file
    // [[codeblock]]
    void compareWithReference(std::vector<Scalar> referenceData)
    {
        // for (const auto dirIdx : directions_)
        // {
        //     const auto K = EffectivePermeability_[dirIdx];
        //     static const Scalar eps = getParam<Scalar>("Problem.TestEpsilon", 1e-3);
        //     if (Dune::FloatCmp::ne<Scalar>(K, referenceData[dirIdx], eps))
        //     {
        //         std::cerr << "Calculated permeability of " << K << " in "
        //                 <<dirName_[dirIdx]<<"-direction does not match with reference value of "
        //                 << referenceData[dirIdx] << std::endl;
        //     }
        // }
    }
    // [[/codeblock]]
    //
    // ### Set the directions that need to be considered
    //
    // [[codeblock]]
    void setDirections(std::vector<std::size_t> directions)
    {
        directions_ = directions;
    }
    // [[/codeblock]]

    // ### Save the relevant data for plot of permeability ratio vs. Forchheimer number
    // [[codeblock]]
private:
    // void writePermeabilityratioVsForchheimerNumber_(std::size_t dirIdx)
    // {
    //     // open a logfile
    //     std::ofstream logfile(dirName_[dirIdx] + "-dir-PermeabilityratioVsForchheimerNumber.dat");

    //     // save the data needed to be plotted in logfile
    //     for (int i = 0; i < apparentPermeability_[dirIdx].size(); i++)
    //     {
    //         // compute the Forchheimer number
    //         const Scalar forchheimerNumber = EffectivePermeability_[dirIdx] * forchheimerCoefficient_[dirIdx] * inertiaToViscousRatio_[dirIdx][i];
    //         // ratio between apparrent permeability and darcy permeability
    //         const Scalar permeabilityRatio = apparentPermeability_[dirIdx][i] / EffectivePermeability_[dirIdx];

    //         logfile << forchheimerNumber << " " << permeabilityRatio << std::endl;
    //     }
    // }
    void writeEffPermeabilityVsAverageSaturation_(std::size_t dirIdx, int phaseIdx)
    {
        // open a logfile
        std::string fileName = dirName_[dirIdx] + "-dir-" + phaseName_[phaseIdx] +"-EffPermeabilityVsAverageSaturation.dat";
        std::ofstream logfile(fileName);

        // save the data needed to be plotted in logfile
        for (int i = 0; i < EffectivePermeability_[dirIdx].size(); i++)
            logfile << averageSaturation_[dirIdx][i] << " " << EffectivePermeability_[dirIdx][i] << std::endl;

        filesToPlot_[0].push_back(fileName);
    }

    void writeCapillaryPressureVsAverageSaturation_(std::size_t dirIdx)
    {
        if(filesToPlot_[1].size())
            return;
        // open a logfile
        std::string fileName = dirName_[dirIdx] + "-dir-CapillaryPressureVsAverageSaturation.dat";

        std::ofstream logfile(fileName);

        // save the data needed to be plotted in logfile
        for (int i = 0; i < EffectivePermeability_[dirIdx].size(); i++)
            logfile << averageSaturation_[dirIdx][i] << " " << pc_[dirIdx][i] << std::endl;

        filesToPlot_[1].push_back(fileName);
    }
    // [[/codeblock]]

    // ### Save the relevant data for plot of inverse of apparent permability vs. rho v / mu
    // [[codeblock]]
    // void writeInversePrmeabilityVsInertiaToViscousRatio_(std::size_t dirIdx)
    // {
    //     // open a logfile and write inverese of apparent permeability given by the model vs. inertial to viscous ratio (rho v / mu)
    //     std::ofstream logfile(dirName_[dirIdx] + "-dir-InversePrmeabilityVsInertiaToViscousRatio.dat");

    //     // save the data needed to be plotted in logfile
    //     for (int i = 0; i < apparentPermeability_[dirIdx].size(); i++)
    //     {
    //         const Scalar inertiaToViscousRatio = inertiaToViscousRatio_[dirIdx][i];
    //         const Scalar inverseAppPermeability = 1 / apparentPermeability_[dirIdx][i];

    //         // compute inverse of apparent permeability using the Forchheimer permeability and coefficient
    //         const Scalar inverseAppPermeabilityForchheimer = 1 / forchheimerPermeability_[dirIdx] + inertiaToViscousRatio * forchheimerCoefficient_[dirIdx];

    //         logfile << inertiaToViscousRatio << " " << 1e-12 * inverseAppPermeability << " " << 1e-12 * inverseAppPermeabilityForchheimer << std::endl;
    //     }
    // }
    // [[/codeblock]]
    //
    // ### Plot permeability ratio vs. Forchheimer number using Gnuplot
    //
    // [[codeblock]]
    static void plotRelPermeabilityVsSw_()
    {
        // using gnuplot interface
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        gnuplot.resetAll();
        std::string title{}, option{};

        int count = 1;
        using std::to_string;
        for (const auto& fileName : filesToPlot_[0])
        {
            std::size_t isWetting = fileName.find("Wetting");

            std::string legend = isWetting != std::string::npos ? "u 1:2 title \" Wetting\" with lines" :
                                                                  "u 1:2 title \" Non-wetting\" with lines";
            // add the data in each direction for plot
            gnuplot.addFileToPlot(fileName, legend);

            // set the properties of lines to be plotted
            option += "set linetype " + to_string(count) + " linecolor " + to_string(count) + " linewidth 5\n";
            count++;
        }

        // report the darcy permeability in each direction as the title of the plot
        title += Fmt::format("Effective permeability vs. Sw");

        option += "set title \"" + title + "\"\n";
        // option += "set logscale x""\n";
        // option += "set format x '10^{%L}'""\n";

        gnuplot.setXlabel("Sw [-]");
        gnuplot.setYlabel("Effective Permeability [m2]");
        gnuplot.setOption(option);
        gnuplot.plot("Effective Permeability vs. Sw");
    }

    static void plotCapillaryPressureVsSw_()
    {
        // using gnuplot interface
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        gnuplot.resetAll();
        std::string title{}, option{};
        for (const auto& fileName : filesToPlot_[1])
        {
            std::string legend = "u 1:2 title \" \" with lines";
            // add the data in each direction for plot
            gnuplot.addFileToPlot(fileName, legend);

            // set the properties of lines to be plotted
            option += "set linetype 1 linecolor 1 linewidth 5\n";
        }

        // report the darcy permeability in each direction as the title of the plot
        title += Fmt::format("Capillary Pressure vs. Sw");

        option += "set title \"" + title + "\"\n";
        // option += "set logscale x""\n";
        // option += "set format x '10^{%L}'""\n";

        gnuplot.setXlabel("Sw [-]");
        gnuplot.setYlabel("Capillary Pressure [Pa]");
        gnuplot.setOption(option);
        gnuplot.plot("Capillary Pressure vs. Sw");
    }
    // [[/codeblock]]
    // [[details]] private data members
    // [[codeblock]]
    std::array<std::vector<Scalar>, 3> samplePointsX_;
    std::array<std::vector<Scalar>, 3> samplePointsY_;
    std::array<std::vector<Scalar>, 3> EffectivePermeability_;
    std::array<std::vector<Scalar>, 3> averageSaturation_;
    std::array<std::vector<Scalar>, 3> pc_;
    const std::array<std::string, 3> dirName_ = {"X", "Y", "Z"};
    const std::array<std::string, 2> phaseName_{"Wetting", "Non-wetting"};
    std::vector<std::size_t> directions_;
    inline static std::array<std::vector<std::string>, 2> filesToPlot_; // 0 index for rel perm vs. sw and 1 index for pc vs. sw
};

} // end namespace Dumux
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
