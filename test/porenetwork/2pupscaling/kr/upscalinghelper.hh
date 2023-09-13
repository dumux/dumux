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
#include <string>

#include <dune/common/float_cmp.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/format.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/common/utilities.hh>
#include <dumux/porenetwork/common/labels.hh>


namespace Dumux::PoreNetwork {

template<class Problem, class GridVariables, class SolutionVector>
class UpscalingHelperTwoP
{
    using Scalar = typename GridVariables::VolumeVariables::PrimaryVariables::value_type;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using FS = typename GridVariables::VolumeVariables::FluidSystem;
    using NumEqVector = typename SolutionVector::block_type;

    struct BoundaryAvgValues
    {
        std::array<Scalar, 2> p;
        std::array<Scalar, 2> density;
        std::array<Scalar, 2> viscosity;
    };
public:

    UpscalingHelperTwoP(const Problem& problem, const GridVariables& gridVariables, const SolutionVector& sol)
    : avgValues_(AveragedValues<GridVariables, SolutionVector>(gridVariables, sol))
    , gridVariables_(gridVariables)
    , sol_(sol)
    , problem_(problem)
    {
        setAvgVariables_();
        logfile_.open("logfile_" + problem_.name() + ".txt"); //for the logfile
        logfile_ <<"Logfile for: " + problem_.name()  << std::endl;
        logfile_ << std::left << std::setw(20) << std::setfill(' ') << "swAveraged"
                 << std::left << std::setw(20) << std::setfill(' ') << "krw"
                 << std::left << std::setw(20) << std::setfill(' ') << "krn"
                 << std::endl;
    }

    void setDataPoints()
    {
        computeEffPerm_();
    }

    template <class LeafGridView, class BoundaryFlux>
    void setDataPoints(const LeafGridView& leafGridView, const BoundaryFlux& boundaryFlux)
    {
        auto inletPoreLabel = problem_.inletPoreLabel();
        auto outletPoreLabel = problem_.outletPoreLabel();

        const auto inletTotalMassFlux = boundaryFlux.getFlux(std::vector<int>{inletPoreLabel}).totalFlux;
        const auto outletTotalMassFlux = boundaryFlux.getFlux(std::vector<int>{outletPoreLabel}).totalFlux;

        computeBoundaryAvgValues_(leafGridView);
        computeVolumeFlux_(inletTotalMassFlux, outletTotalMassFlux);
        computeEffPerm_();
    }

    template <class LeafGridView, class BoundaryFlux>
    bool reachedEquilibrium(const LeafGridView& leafGridView, const BoundaryFlux& boundaryFlux)
    {
        auto inletPoreLabel = problem_.inletPoreLabel();
        auto outletPoreLabel = problem_.outletPoreLabel();

        const auto inletTotalMassFlux = boundaryFlux.getFlux(std::vector<int>{inletPoreLabel}).totalFlux;
        const auto outletTotalMassFlux = boundaryFlux.getFlux(std::vector<int>{outletPoreLabel}).totalFlux;

        computeBoundaryAvgValues_(leafGridView);
        computeVolumeFlux_(inletTotalMassFlux, outletTotalMassFlux);
        return checkIfInEquilibrium_();
    }

    template <class LeafGridView>
    bool reachedEquilibrium(const LeafGridView& leafGridView, Scalar dt)
    {

        computeDomainAvgValues_(leafGridView);
        return checkIfInEquilibrium_(dt);
    }

    // ### Set data points to calculate intrinsic permeability and Forchheimer coefficient
    // This function first evaluates the mass flux leaving the network in the direction of the applied pressure gradient.
    // Afterwards, the mass flux is converted into an volume flux which is used to calculate the apparent velocity.
    // Then apparent permeability of the network is computed and stored for furthure calculations.
    // [[codeblock]]
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
    // [[codeblock]]


    std::array<Scalar, 2> volumeFluxes(const std::string& boundaryName) const
    {
        if (boundaryName == "inlet")
            return inletVolumeFlux_;

        return outletVolumeFlux_;
    }

    const Scalar effectivePermeability(int phaseIdx) const
    {   return effPerm_[phaseIdx];  }

    bool checkIfInEquilibrium()
    {   return checkIfInEquilibrium_(); }

    // ### Determine the domain's side lengths
    //
    // We determine the domain side length by using the bounding box of the network
    // [[codeblock]]
    template<class GridGeometry>
    auto setSideLengths(const GridGeometry& gridGeometry)
    {
        std::cout << "Automatically determining side lengths of REV based on bounding box of pore network" << std::endl;
        for (int dimIdx = 0; dimIdx < GridGeometry::GridView::dimensionworld; ++dimIdx)
            length_[dimIdx] = gridGeometry.bBoxMax()[dimIdx] - gridGeometry.bBoxMin()[dimIdx];

        return length_;
    }
    // [[/codeblock]]


private:

    void setAvgVariables_()
    {
        avgValues_.addAveragedQuantity([](const auto& v){ return v.saturation(FS::phase0Idx); }, [](const auto& v){ return v.poreVolume(); }, "avgSat");

        avgValues_.addAveragedQuantity([](const auto& v){ return v.pressure(FS::phase0Idx); }, [](const auto& v){ return v.saturation(FS::phase0Idx)*v.poreVolume(); }, "avgPw");
        avgValues_.addAveragedQuantity([](const auto& v){ return v.pressure(FS::phase1Idx); }, [](const auto& v){ return v.saturation(FS::phase1Idx)*v.poreVolume(); }, "avgPn");

        avgValues_.addAveragedQuantity([](const auto& v){ return v.density(FS::phase0Idx); }, [](const auto& v){ return v.saturation(FS::phase0Idx)*v.poreVolume(); }, "avgDensityw");
        avgValues_.addAveragedQuantity([](const auto& v){ return v.density(FS::phase1Idx); }, [](const auto& v){ return v.saturation(FS::phase1Idx)*v.poreVolume(); }, "avgDensityn");

        avgValues_.addAveragedQuantity([](const auto& v){ return v.viscosity(FS::phase0Idx); }, [](const auto& v){ return v.saturation(FS::phase0Idx)*v.poreVolume(); }, "avgViscosityw");
        avgValues_.addAveragedQuantity([](const auto& v){ return v.viscosity(FS::phase1Idx); }, [](const auto& v){ return v.saturation(FS::phase1Idx)*v.poreVolume(); }, "avgViscosityn");
    }

    template <class LeafGridView>
    void  computeBoundaryAvgValues_(const LeafGridView& leafGridView)
    {
        auto dofsToNeglect = dofsToNeglect_(leafGridView, std::vector<int>{Labels::outlet, Labels::interior});

        avgValues_.eval(dofsToNeglect);
        assignBoundaryAvgValues_(inletAvgValues_, avgValues_);

        dofsToNeglect.clear();
        dofsToNeglect = dofsToNeglect_(leafGridView, std::vector<int>{Labels::inlet, Labels::interior});

        avgValues_.eval(dofsToNeglect);
        assignBoundaryAvgValues_(outletAvgValues_, avgValues_);
    }

    template <class LeafGridView>
    void computeDomainAvgValues_(const LeafGridView& leafGridView)
    {
        auto dofsToNeglect = dofsToNeglect_(leafGridView, std::vector<int>{Labels::outlet/*, Labels::inlet*/});

        avgValues_.eval(dofsToNeglect);
        const Scalar avgSw = avgValues_["avgSat"];

        // store the three most recent averaged saturations
        std::rotate(swAvg_.rbegin(), swAvg_.rbegin()+1, swAvg_.rend());
        swAvg_[0]= avgSw;
    }

    void assignBoundaryAvgValues_(BoundaryAvgValues& boundaryAvgValues, const AveragedValues<GridVariables, SolutionVector>& avgValues)
    {
        boundaryAvgValues.p[FS::phase0Idx] = avgValues["avgPw"];
        boundaryAvgValues.p[FS::phase1Idx] = avgValues["avgPn"];
        boundaryAvgValues.density[FS::phase0Idx] = avgValues["avgDensityw"];
        boundaryAvgValues.density[FS::phase1Idx] = avgValues["avgDensityn"];
        boundaryAvgValues.viscosity[FS::phase0Idx] = avgValues["avgViscosityw"];
        boundaryAvgValues.viscosity[FS::phase1Idx] = avgValues["avgViscosityn"];
    }

    template <class LeafGridView>
    std::vector<std::size_t> dofsToNeglect_(const LeafGridView& leafGridView, std::vector<int> labels)
    {
        std::vector<std::size_t> dofsToNeglect;
        for (const auto& vertex : vertices(leafGridView))
        {
            const auto vIdx = problem_.gridGeometry().vertexMapper().index(vertex);
            const auto poreLabel = problem_.gridGeometry().poreLabel(vIdx);
            if (std::any_of(labels.begin(), labels.end(), [&](const int l){ return l == poreLabel; }))
                dofsToNeglect.push_back(vIdx);
        }
        return dofsToNeglect;
    }

    void computeVolumeFlux_(const NumEqVector& inletTotalMassFlux, const NumEqVector& outletTotalMassFlux)
    {
        inletVolumeFlux_[FS::phase0Idx] = inletTotalMassFlux[FS::phase0Idx]/inletAvgValues_.density[FS::phase0Idx];
        inletVolumeFlux_[FS::phase1Idx] = inletTotalMassFlux[FS::phase1Idx]/inletAvgValues_.density[FS::phase1Idx];
        outletVolumeFlux_[FS::phase0Idx] = outletTotalMassFlux[FS::phase0Idx]/outletAvgValues_.density[FS::phase0Idx];
        outletVolumeFlux_[FS::phase1Idx] = outletTotalMassFlux[FS::phase1Idx]/outletAvgValues_.density[FS::phase1Idx];
    }

    void computeEffPerm_()
    {

        // get the domain side lengths from the problem
        auto sideLengths = length_;
        auto currentDirection = problem_.direction();
        Scalar L = sideLengths[currentDirection];
        // calculate apparent velocity
        sideLengths[currentDirection] = 1.0;
        const auto outflowArea = std::accumulate(sideLengths.begin(), sideLengths.end(), 1.0, std::multiplies<Scalar>());
        for (int phaseIdx = 0; phaseIdx < 2; phaseIdx++)
        {
            const Scalar mu = outletAvgValues_.viscosity[phaseIdx];
            effPerm_[currentDirection][phaseIdx] = outletVolumeFlux_[phaseIdx] * mu * L /(outflowArea * (inletAvgValues_.p[phaseIdx] - outletAvgValues_.p[phaseIdx]));
        }

        logfile_ << std::scientific << std::setprecision(2)<< std::left << std::setw(20) << std::setfill(' ') << swAvg_[0]
                         << std::left << std::setw(20) << std::setfill(' ') <<  effPerm_[currentDirection][0]
                         << std::left << std::setw(20) << std::setfill(' ') << effPerm_[currentDirection][1]
                         << std::endl;
        std::cout<<"   effPerm_[direction_]   "<<effPerm_[currentDirection][0]<<std::endl;
    }

    bool checkIfInEquilibrium_()
    {
        using std::abs;
        Scalar diff = 0.0;
        for (int i = 0; i<2 ; i++)
        {
            diff = (abs(inletVolumeFlux_[i]) - abs(outletVolumeFlux_[i]));
            if (abs(diff) > abs(outletVolumeFlux_[i]) * 1e-2)
                return 0;
        }

        return 1;
    }

    bool checkIfInEquilibrium_(Scalar dt)
    {
        using std::abs;
        Scalar dSwDt = abs(swAvg_[0]-swAvg_[1]) / dt;
        std::cout<<"  equilibrium????   "<< dSwDt<<"   "<<swAvg_[0]<<"   "<<swAvg_[1]<<std::endl;
        if (dSwDt < 1e-15)
        {
            std::cout<<"  equilibrium????   "<< dSwDt<<"   "<<swAvg_[0]<<"   "<<swAvg_[1]<<std::endl;
            swAvg_[1] = 2.0;
            return 1;
        }

        return 0;
    }



    // ### Save the relevant data for plot of permeability ratio vs. Forchheimer number
    // [[codeblock]]
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
    AveragedValues<GridVariables, SolutionVector> avgValues_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    const Problem& problem_;
    BoundaryAvgValues inletAvgValues_;
    BoundaryAvgValues outletAvgValues_;
    std::array<Scalar, 2> inletVolumeFlux_;
    std::array<Scalar, 2> outletVolumeFlux_;
    std::array<std::array<Scalar, 2>, 3>  effPerm_;
    std::array<Scalar, 2> swAvg_ = {{1.0, 1.0}};
    std::array<Scalar, 3> length_;
    std::ofstream logfile_;

};

} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
