// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
#ifndef DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH
#define DUMUX_PLOT_THERMAL_CONDUCTIVITY_LAW_HH

#include <string>
#include <vector>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux {

// forward declaration
template<class Scalar> class GnuplotInterface;

/*!
 * \ingroup InputOutput
 * \brief Interface for plotting the non-isothermal two-phase fluid-matrix-interaction laws
 */
template<class Scalar, class ThermalConductivityModel, class FS>
class PlotThermalConductivityModel
{
    using FluidSystem = FS;
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;

    // phase indices
    enum
    {
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx
    };

public:
    /*!
     * \brief Constructor
     *
     * Initializes the fluid system.
     *
     * \param temperature temperature in \f$\mathrm{[K]}\f$
     * \param pressure reference pressure in \f$\mathrm{[Pa]}\f$
     */
    PlotThermalConductivityModel(Scalar temperature = 283.15,
                                 Scalar pressure = 1e5)
    : numIntervals_(1000)
    {
        FluidState fluidstate;
        fluidstate.setTemperature(temperature);
        fluidstate.setPressure(phase0Idx, pressure);
        fluidstate.setPressure(phase1Idx, pressure);
        lambdaW_ = FluidSystem::thermalConductivity(fluidstate, phase0Idx);
        lambdaN_ = FluidSystem::thermalConductivity(fluidstate, phase1Idx);
    }

    /*!
     * \brief Add a effective thermal conductivity-saturation curve to the plot
     *
     * \param gnuplot The gnuplot interface
     * \param porosity The porosity
     * \param rhoSolid The density of the solid phase
     * \param lambdaSolid The conductivity of the solid phase
     * \param lowerSat Minimum x-value for data set
     * \param upperSat Maximum x-value for data set
     * \param curveName Name of the data set
     * \param curveOptions Plotting options associated with that data set
     */
    void addlambdaeffcurve(GnuplotInterface<Scalar> &gnuplot,
                           Scalar porosity,
                           Scalar rhoSolid,
                           Scalar lambdaSolid,
                           Scalar lowerSat = 0.0,
                           Scalar upperSat = 1.0,
                           std::string curveName = "lambdaeff",
                           std::string curveOptions = "w l")
    {
        std::vector<Scalar> sw(numIntervals_+1);
        std::vector<Scalar> lambda(numIntervals_+1);
        Scalar satInterval = upperSat - lowerSat;

        for (int i = 0; i <= numIntervals_; i++)
        {
            sw[i] = lowerSat + satInterval * Scalar(i) / Scalar(numIntervals_);
            VolumeVariables volVars(sw[i], lambdaN_, lambdaW_, lambdaSolid, porosity, rhoSolid);
            lambda[i] = ThermalConductivityModel::effectiveThermalConductivity(volVars);
        }

        gnuplot.addDataSetToPlot(sw, lambda, curveName, curveOptions);
    }

private:

    class VolumeVariables
    {
    public:
        VolumeVariables(Scalar saturation, Scalar lambdaN, Scalar lambdaW, Scalar lambdaSolid, Scalar porosity, Scalar rhoSolid)
        : saturation_(saturation)
        , lambdaN_(lambdaN)
        , lambdaW_(lambdaW)
        , lambdaSolid_(lambdaSolid)
        , porosity_(porosity)
        , rhoSolid_(rhoSolid)
        {}

        using FluidSystem = typename PlotThermalConductivityModel::FluidSystem;

        Scalar saturation(const int phaseIdx) const
        {
            if (phaseIdx == wettingPhase())
                return saturation_;
            else
                return 1.0 - saturation_;
        }

        Scalar fluidThermalConductivity(const int phaseIdx) const
        {
            if (phaseIdx == wettingPhase())
                return lambdaW_;
            else
                return lambdaN_;
        }

        int wettingPhase() const
        { return phase0Idx; }

        Scalar porosity() const
        { return porosity_; }

        Scalar solidThermalConductivity() const
        { return lambdaSolid_; }

        Scalar solidDensity() const
        { return rhoSolid_;}

    private:
        Scalar saturation_;
        Scalar lambdaN_;
        Scalar lambdaW_;
        Scalar lambdaSolid_;
        Scalar porosity_;
        Scalar rhoSolid_;
    };

    int numIntervals_;
    Scalar lambdaN_;
    Scalar lambdaW_;
};

} // end namespace Dumux

#endif
