// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/format.hh>

namespace Dumux {

template<class Scalar>
class UpscalingHelper
{
public:

    // ### Set sample points to calculate intrinsic permeability and Forchheimer coefficient
    // This function first evaluates the mass flux leaving the network in the direction of the applied pressure gradient.
    // Afterwards, the mass flux is converted into an area specify volume flux which with its corresponding pressure gradient are stored as
    // sample points to be used in regression operation to find intrinsic permeability and Forchheimer coefficient
    // [[codeblock]]
    template<class Problem>
    void setSamplePoints(const Problem& problem, const Scalar totalMassFlux)
    {
        // get the domain side lengths from the problem
        auto sideLengths = problem.sideLengths();

        // convert mass to volume flux
        const auto volumeFlux = totalMassFlux / problem.liquidDensity();;

        // calculate apparent velocity
        sideLengths[problem.direction()] = 1.0;
        const auto outflowArea = std::accumulate(sideLengths.begin(), sideLengths.end(), 1.0, std::multiplies<Scalar>());
        const auto vApparent= volumeFlux / outflowArea;

        // set sample point for permability calculation
        const auto samplePointY = problem.pressureGradient() / problem.liquidDynamicViscosity() / vApparent;
        const auto samplePointX = problem.liquidDensity() * vApparent / problem.liquidDynamicViscosity();

        samplePointsX_[problem.direction()].push_back(samplePointX);
        samplePointsY_[problem.direction()].push_back(samplePointY);

        // compute apparent permeability
        const auto K = vApparent / problem.pressureGradient() * problem.liquidDynamicViscosity();

        // calculate Forchheimer number (Forchheimer coefficient will be included later)
        const auto forchheimerNumber = problem.liquidDensity() * vApparent / problem.liquidDynamicViscosity();

        // store apparent permeability and corresponding Forchheimer number
        apparentPermeability_[problem.direction()].push_back(K);
        forchheimerNumber_[problem.direction()].push_back(forchheimerNumber);
    }
    // [[/codeblock]]

    // ### Calculate intrinsic permeability and Forchheimer coefficient.
    // This function first calculate intrinsic permeability and Forchheimer coefficient using linear least squares regression method
    // and reports them. It also plot the apparent permeability of the porous medium versus Forchheimer number/pressure gradient in each
    // simulation.
    // [[codeblock]]
    void calculateUpscaledProperties(bool isCreepingFlow)
    {
        for (int dirIdx = 0; dirIdx < 3; dirIdx++)
        {
            // determine Darcy permeability as the maximum permeability of the domain
            darcyPermeability_[dirIdx] = *max_element(apparentPermeability_[dirIdx].begin(), apparentPermeability_[dirIdx].end());
            if (!isCreepingFlow)
            {
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
        // using gnuplot interface
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        std::string title{}, option{};
        for (int dirIdx = 0; dirIdx < 3; dirIdx++)
        {
            // add the data in each direction for plot
            gnuplot.addFileToPlot(dirName_[dirIdx] + "-dir.dat");
            // set the properties of lines to be plotted
            option += Fmt::format("set linetype {0} linecolor {0} linewidth 7\n", dirIdx+1);
            // report the darcy permeability in each direction as the title of the plot
            title += Fmt::format("{}-permeability= {:.3e} m^2   ", dirName_[dirIdx], darcyPermeability_[dirIdx]);
        }
        option +="set title \"" + title + "\"\n";
        option += "set logscale x""\n";
        option += "set format x '10^{%L}'""\n";

        gnuplot.setXlabel("Forchheimer Number [-]");
        gnuplot.setYlabel("Apparent permeability / Darcy permeability [-]");
        gnuplot.setOption(option);
        gnuplot.plot("permeability_ratio_versus_forchheimer_number");
    }
    // [[/codeblock]]
    //
    // ### Save the relevant data for plot
    // [[codeblock]]
    void writePlotDataToFile(std::size_t dirIdx)
    {
        // Open a logfile
        std::ofstream logfile(dirName_[dirIdx]+"-dir.dat");

        // Save the data needed to be plotted in logfile
        for (int i = 0; i < apparentPermeability_[dirIdx].size(); i++)
        {
             // Include characteristics length, sqrt(permeability) to Reynolds number calculation
            const Scalar forchheimerNumber
                = darcyPermeability_[dirIdx] * forchheimerCoefficient_[dirIdx] * forchheimerNumber_[dirIdx][i];
            // Ratio between apparrent permeability and darcy permeability
            const Scalar permeabilityRatio
                = apparentPermeability_[dirIdx][i] / darcyPermeability_[dirIdx];

            logfile << forchheimerNumber<< " " << permeabilityRatio << std::endl;
        }
    }
    // [[/codeblock]]
    //
    // ### Report the upscaled data
    // [[codeblock]]
    void report(bool isCreepingFlow)
    {
        // Report the results for each direction
        for (int dirIdx = 0; dirIdx < 3; dirIdx++)
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
        for (int dirIdx = 0; dirIdx < 3; dirIdx++)
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
private:
    std::array<std::vector<Scalar>, 3> samplePointsX_;
    std::array<std::vector<Scalar>, 3> samplePointsY_;
    std::array<std::vector<Scalar>, 3> apparentPermeability_;
    std::array<std::vector<Scalar>, 3> forchheimerNumber_;
    std::array<Scalar, 3> darcyPermeability_;
    std::array<Scalar, 3> forchheimerPermeability_;
    std::array<Scalar, 3> forchheimerCoefficient_;
    const std::array<std::string, 3> dirName_ = {"X", "Y", "Z"};
};

} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
