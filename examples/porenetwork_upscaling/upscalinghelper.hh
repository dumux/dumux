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
    template <class Problem>
    void setDataPoints(const Problem &problem, const Scalar totalMassFlux)
    {
        // get the domain side lengths from the problem
        auto sideLengths = problem.sideLengths();


        // get the applied pressure gradient
        const auto pressureGradient = problem.pressureGradient();
        const auto pressureDrop = pressureGradient * sideLengths[problem.direction()];

        // get the fluid properties
        const auto liquidDensity = problem.liquidDensity();
        const auto liquidDynamicViscosity = problem.liquidDynamicViscosity();

        // convert mass to volume flux
        const auto volumeFlux = totalMassFlux / liquidDensity;
        ;

        // calculate apparent velocity
        sideLengths[problem.direction()] = 1.0;
        const auto outflowArea = std::accumulate(sideLengths.begin(), sideLengths.end(), 1.0, std::multiplies<Scalar>());
        const auto vApparent = volumeFlux / outflowArea;

        // compute apparent permeability
        const auto KApparent = vApparent / pressureGradient * liquidDynamicViscosity;
        // calculate rho v / mu, called inertia to viscous ratio in the rest of the code
        const auto inertiaToViscousRatio = liquidDensity * vApparent / liquidDynamicViscosity;

        // store the required data for further calculations
        totalPressureDrop_[problem.direction()].push_back(pressureDrop);
        apparentVelocity_[problem.direction()].push_back(vApparent);
        apparentPermeability_[problem.direction()].push_back(KApparent);
        inertiaToViscousRatio_[problem.direction()].push_back(inertiaToViscousRatio);
    }
    // [[/codeblock]]

    // ### Calculate intrinsic permeability and Forchheimer coefficient.
    // This function first calculate intrinsic permeability and Forchheimer coefficient using linear least squares regression method
    // and reports them. It also plot the apparent permeability of the porous medium versus Forchheimer number/pressure gradient in each
    // simulation.
    // [[codeblock]]
    template <class Problem>
    void calculateUpscaledProperties(const Problem &problem, bool isCreepingFlow)
    {
        const auto sideLengths = problem.sideLengths();
        const auto liquidDynamicViscosity = problem.liquidDynamicViscosity();

        for (const auto dirIdx : directions_)
        {
            // determine Darcy permeability as the maximum permeability of the domain
            darcyPermeability_[dirIdx] = *max_element(apparentPermeability_[dirIdx].begin(), apparentPermeability_[dirIdx].end());
            if (!isCreepingFlow)
            {
                for (int i = 0; i < totalPressureDrop_[dirIdx].size(); i++)
                {
                    // calculate the Darcy pressure drop.
                    const Scalar darcyPressureDrop = liquidDynamicViscosity * apparentVelocity_[dirIdx][i] * sideLengths[dirIdx] / darcyPermeability_[dirIdx];

                    // calculate the ratio of Dracy to total pressure drop
                    const Scalar pressureDropRatio = darcyPressureDrop / totalPressureDrop_[dirIdx][i];

                    // set sample points for upscaling of Forchheimer parameters.
                    // first, check the permability ratio to see if the flow regime is Forchheimer.
                    if (pressureDropRatio < 0.99)
                    {
                        samplePointsX_[dirIdx].push_back(inertiaToViscousRatio_[dirIdx][i]);
                        samplePointsY_[dirIdx].push_back(1 / apparentPermeability_[dirIdx][i]);
                    }
                }
                // determine regression line and accordingly the Forchheimer permeability and the Forchheimer coefficient
                const auto [intercept, slope] = linearRegression(samplePointsX_[dirIdx], samplePointsY_[dirIdx]);
                forchheimerPermeability_[dirIdx] = 1.0 / intercept;
                forchheimerCoefficient_[dirIdx] = slope;
                writePlotDataToFile(dirIdx);
            }
        }
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
    void plot()
    {
        // plot permeability ratio vs. Forchheimer number
        plotPermeabilityratioVsForchheimerNumber_();

        // plot inverse of apparent permability vs. rho v / mu
        plotInversePrmeabilityVsInertiaToViscousRatio_();
    }
    // [[/codeblock]]
    //
    // ### Save the relevant data for plot
    // [[codeblock]]
    void writePlotDataToFile(std::size_t dirIdx)
    {
        // write permeability ratio vs. Forchheimer number
        writePermeabilityratioVsForchheimerNumber_(dirIdx);

        // write inverse of apparent permability vs. rho v / mu
        writeInversePrmeabilityVsInertiaToViscousRatio_(dirIdx);
    }
    // [[/codeblock]]
    //
    // ### Report the upscaled data
    // [[codeblock]]
    void report(bool isCreepingFlow)
    {
        // Report the results for each direction
        for (const auto dirIdx : directions_)
        {
            std::cout << Fmt::format("\n{:#>{}}\n\n", "", 40)
                      << Fmt::format("{}-direction:\n", dirName_[dirIdx])
                      << Fmt::format("-- Darcy (intrinsic) permeability = {:.3e} m^2\n", darcyPermeability_[dirIdx]);

            // Report non-creeping flow upscaled properties
            if (!isCreepingFlow)
            {
                std::cout << Fmt::format("-- Forchheimer permeability = {:.3e} m^2\n", forchheimerPermeability_[dirIdx]);
                std::cout << Fmt::format("-- Forchheimer coefficient = {:.3e} m^-1\n", forchheimerCoefficient_[dirIdx]);
            }

            std::cout << Fmt::format("\n{:#>{}}\n", "", 40) << std::endl;
        }
    }
    // [[/codeblock]]
    //
    // ### Compare with reference data provided in input file
    // [[codeblock]]
    void compareWithReference(std::vector<Scalar> referenceData)
    {
        for (const auto dirIdx : directions_)
        {
            const auto K = darcyPermeability_[dirIdx];
            static const Scalar eps = getParam<Scalar>("Problem.TestEpsilon", 1e-3);
            if (Dune::FloatCmp::ne<Scalar>(K, referenceData[dirIdx], eps))
            {
                std::cerr << "Calculated permeability of " << K << " in "
                        <<dirName_[dirIdx]<<"-direction does not match with reference value of "
                        << referenceData[dirIdx] << std::endl;
            }
        }
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
    void writePermeabilityratioVsForchheimerNumber_(std::size_t dirIdx)
    {
        // open a logfile
        std::ofstream logfile(dirName_[dirIdx] + "-dir-PermeabilityratioVsForchheimerNumber.dat");

        // save the data needed to be plotted in logfile
        for (int i = 0; i < apparentPermeability_[dirIdx].size(); i++)
        {
            // compute the Forchheimer number
            const Scalar forchheimerNumber = darcyPermeability_[dirIdx] * forchheimerCoefficient_[dirIdx] * inertiaToViscousRatio_[dirIdx][i];
            // ratio between apparrent permeability and darcy permeability
            const Scalar permeabilityRatio = apparentPermeability_[dirIdx][i] / darcyPermeability_[dirIdx];

            logfile << forchheimerNumber << " " << permeabilityRatio << std::endl;
        }
    }
    // [[/codeblock]]

    // ### Save the relevant data for plot of inverse of apparent permability vs. rho v / mu
    // [[codeblock]]
    void writeInversePrmeabilityVsInertiaToViscousRatio_(std::size_t dirIdx)
    {
        // open a logfile and write inverese of apparent permeability given by the model vs. inertial to viscous ratio (rho v / mu)
        std::ofstream logfile(dirName_[dirIdx] + "-dir-InversePrmeabilityVsInertiaToViscousRatio.dat");

        // save the data needed to be plotted in logfile
        for (int i = 0; i < apparentPermeability_[dirIdx].size(); i++)
        {
            const Scalar inertiaToViscousRatio = inertiaToViscousRatio_[dirIdx][i];
            const Scalar inverseAppPermeability = 1 / apparentPermeability_[dirIdx][i];

            // compute inverse of apparent permeability using the Forchheimer permeability and coefficient
            const Scalar inverseAppPermeabilityForchheimer = 1 / forchheimerPermeability_[dirIdx] + inertiaToViscousRatio * forchheimerCoefficient_[dirIdx];

            logfile << inertiaToViscousRatio << " " << 1e-12 * inverseAppPermeability << " " << 1e-12 * inverseAppPermeabilityForchheimer << std::endl;
        }
    }
    // [[/codeblock]]
    //
    // ### Plot permeability ratio vs. Forchheimer number using Gnuplot
    //
    // [[codeblock]]
    void plotPermeabilityratioVsForchheimerNumber_()
    {
        // using gnuplot interface
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);

        for (const auto dirIdx : directions_)
        {
            gnuplot.resetAll();
            std::string title{}, option{};

            // add the data in each direction for plot
            gnuplot.addFileToPlot(dirName_[dirIdx] + "-dir-PermeabilityratioVsForchheimerNumber.dat", "notitle with lines");
            // set the properties of lines to be plotted
            option += "set linetype 1 linecolor 1 linewidth 7\n";
            // report the darcy permeability in each direction as the title of the plot
            title += Fmt::format("{}-direction, Darcy permeability= {:.3e} m^2   ", dirName_[dirIdx], darcyPermeability_[dirIdx]);

            option += "set title \"" + title + "\"\n";
            option += "set logscale x""\n";
            option += "set format x '10^{%L}'""\n";

            gnuplot.setXlabel("Forchheimer Number [-]");
            gnuplot.setYlabel("Apparent permeability / Darcy permeability [-]");
            gnuplot.setOption(option);
            gnuplot.plot("permeability_ratio_versus_forchheimer_number");
        }
    }
    // [[/codeblock]]
    //
    // ### Plot inverse of apparent permability vs. rho v / mu using Gnuplot
    //
    // [[codeblock]]
    void plotInversePrmeabilityVsInertiaToViscousRatio_()
    {
        // using gnuplot interface
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);

        for (const auto dirIdx : directions_)
        {
            gnuplot.resetAll();
            std::string title{}, option{};
            std::string legend0 = "u 1:2 title \"Network model\" with lines";
            // add the data in each direction for plot, first set of data
            gnuplot.addFileToPlot(dirName_[dirIdx] + "-dir-InversePrmeabilityVsInertiaToViscousRatio.dat", legend0);

            // set the properties of lines to be plotted
            option += "set linetype 1 linecolor 1 linewidth 5\n";

            std::string legend1 = "u 1:3 title \"Forchheimer equation\" with lines";
            // add the data in each direction for plot, second set of data
            gnuplot.addFileToPlot(dirName_[dirIdx] + "-dir-InversePrmeabilityVsInertiaToViscousRatio.dat", legend1);

            // set the properties of lines to be plotted
            option += "set linetype 2 linecolor 2 linewidth 5\n";

            // report the darcy permeability in each direction as the title of the plot
            title += Fmt::format("{}-direction, Darcy permeability= {:.3e} m^2   ", dirName_[dirIdx], darcyPermeability_[dirIdx]);

            option += "set title \"" + title + "\"\n";
            option += "set logscale x""\n";
            option += "set format x '10^{%L}'""\n";

            gnuplot.setXlabel("{/Symbol r} v / {/Symbol m} [1/m]");
            gnuplot.setYlabel("1/ Apparent permeability [1/m^2]  x 1e12");
            gnuplot.setOption(option);
            gnuplot.plot("inverse_apppermeability_versus_rhovmu-1");
        }
    }
    // [[/codeblock]]
    // [[details]] private data members
    // [[codeblock]]
    std::array<std::vector<Scalar>, 3> samplePointsX_;
    std::array<std::vector<Scalar>, 3> samplePointsY_;
    std::array<std::vector<Scalar>, 3>totalPressureDrop_;
    std::array<std::vector<Scalar>, 3> apparentVelocity_;
    std::array<std::vector<Scalar>, 3> apparentPermeability_;
    std::array<std::vector<Scalar>, 3> inertiaToViscousRatio_;
    std::array<Scalar, 3> darcyPermeability_;
    std::array<Scalar, 3> forchheimerPermeability_;
    std::array<Scalar, 3> forchheimerCoefficient_;
    const std::array<std::string, 3> dirName_ = {"X", "Y", "Z"};
    std::vector<std::size_t> directions_;
};

} // end namespace Dumux
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
