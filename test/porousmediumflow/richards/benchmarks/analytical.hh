// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Helper functions for analytical solutions
 */

#ifndef DUMUX_RICHARDS_BENCHMARKS_ANALYTICAL_HH
#define DUMUX_RICHARDS_BENCHMARKS_ANALYTICAL_HH

#include <cmath>
#include <algorithm>

#include <dumux/io/format.hh>
#include <dumux/io/gnuplotinterface.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>

#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/components/simpleh2o.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup RichardsTests
 * \brief Analytical solution for infiltration scenario
 * \note Root-soil benchmark paper Schnepf et al. (case M2.1, Eq. 4) https://doi.org/10.3389/fpls.2020.00316
 * \note based on Vanderborght 2005 (see Fig. 4abc and Eq. 56-60) https://doi.org/10.2113/4.1.206
 *
 * Assumption is that the infiltration front is sufficiently far
 * away from the boundary. The boundary conditions for the analytical
 * solution are prescribed at z = ±∞.
 */
class AnalyticalSolutionM21
{
    using PcKrSwCurve = FluidMatrix::VanGenuchtenNoReg<double>;
public:
    AnalyticalSolutionM21()
    : pcKrSwCurve_("SpatialParams")
    {
        gravity_ = 9.81;
        density_ = Components::SimpleH2O<double>::liquidDensity(0,0);
        viscosity_ = Components::SimpleH2O<double>::liquidViscosity(0,0);
        porosity_ = getParam<double>("SpatialParams.Porosity");
        permeability_ = getParam<double>("SpatialParams.Permeability");

        const auto initialHead = getParam<double>("Problem.InitialHeadInCm", -400.0)*0.01; // cm to m
        const auto pwInitial = initialHead*gravity_*density_ + 1.0e5;
        const auto pcInitial = 1.0e5 - pwInitial;
        const auto swIntial = pcKrSwCurve_.sw(pcInitial);

        // compute θ_i, K_i, θ_sur, K_sur, θ_a
        waterContentInitial_ = swIntial*porosity_;
        conductivityInitial_ = conductivity(waterContentInitial_);
        waterContentSurface_ = getParam<double>("Analytical.WaterContentSurface");
        conductivitySurface_ = conductivity(waterContentSurface_);
        waterContentRef_ = 0.5*(waterContentSurface_ + waterContentInitial_);

        std::cout << std::endl << "Computed analytical solution (Vanderborght 2005 4abc, Schnepf 2020 M2.1) with:\n";
        const auto soilType = getParam<std::string>("SpatialParams.Name");
        std::cout << "Material: " << soilType << "\n";
        std::cout << Fmt::format("-> θ_i: {}\n-> θ_sur: {}\n", waterContentInitial_, waterContentSurface_);
        std::cout << Fmt::format("-> θ_a: {}\n", waterContentRef_);
        std::cout << Fmt::format("-> K_i: {}\n-> K_sur: {}\n", conductivityInitial_, conductivitySurface_);

        // create a table for fast evaluation of the inverse
        // for the inverse we need to go from high water content to low water content so
        // that deltaEta is monotonically increasing which is a requirement for the value key
        // for the table interpolation
        waterContents_ = linspace(waterContentSurface_ - 1e-8, waterContentInitial_ + 1e-8, 10000);
        deltaEtas_.resize(waterContents_.size());
        std::transform(waterContents_.begin(), waterContents_.end(), deltaEtas_.begin(),
                       [this](auto theta){ return deltaEtaExact(theta); });

        // plot delta eta over water content
        GnuplotInterface<double> gnuplot(false);
        gnuplot.setOpenPlotWindow(false);
        gnuplot.resetPlot();
        gnuplot.setXlabel("water content");
        gnuplot.setYlabel("delta eta");
        gnuplot.addDataSetToPlot(waterContents_, deltaEtas_, "theta_vs_deltaeta_" + soilType + ".dat", "w l t 'analytical'");
        gnuplot.plot("infiltration_theta_vs_deltaeta");

        // TODO how to choose these? Is there an exact value?
        waterContentRefDepth_ = getParam<double>("Analytical.WaterContentRefDepth");
        waterContentRefTime_ = getParam<double>("Analytical.WaterContentRefTime");
        std::cout << Fmt::format("-> η_a: {}\n", eta(waterContentRefDepth_, waterContentRefTime_));
    }

    // analytical solution for given x and t
    double waterContent(double x /*in cm*/, double t /*in days*/) const
    {
        const auto deltaEta = eta(x, t) - eta(waterContentRefDepth_, waterContentRefTime_);
        return interpolate<InterpolationPolicy::LinearTable>(deltaEta, deltaEtas_, waterContents_);
    }

    double deltaEta(double x /*in cm*/, double t /*in days*/) const
    { return eta(x, t) - eta(waterContentRefDepth_, waterContentRefTime_); }

    double initialWaterContent() const
    { return waterContentInitial_; }

    double surfaceWaterContent() const
    { return waterContentSurface_; }

    double refWaterContent() const
    { return waterContentRef_; }

    auto waterContentAndDeltaEtaPlot() const
    { return std::tie(waterContents_, deltaEtas_); }

private:
    // hydraulic conductivity Kf (including relative permeability) in cm/day
    double conductivity(double waterContent) const
    {
        const auto sw = waterContent/porosity_;
        const auto krw = pcKrSwCurve_.krw(sw);
        return krw*permeability_*gravity_*density_/viscosity_*100*86400;
    }

    // D = Kf*dh/dtheta Van Genuchten 1980 Eq. 10 in cm^2/day
    double diffusivity(double waterContent) const
    {
        // pc = -h*gravity_*density_
        // h = -pc/(gravity_*density_)
        const auto sw = waterContent/porosity_;
        const auto dh_dpc = -1.0/(0.01*gravity_*density_);
        const auto dpc_dsw = pcKrSwCurve_.dpc_dsw(sw);
        const auto dsw_dtheta = 1.0/porosity_;
        const auto dh_dtheta = dh_dpc*dpc_dsw*dsw_dtheta;
        return conductivity(waterContent)*dh_dtheta;
    }

    // integral travelling wave (Eq. 4 Schnepf 2020)
    double integral(double waterContent) const
    {
        const auto integrand = [this](auto theta){
            return diffusivity(theta)/(
                  (conductivitySurface_ - conductivityInitial_)*(theta - waterContentInitial_)
                - (conductivity(theta) - conductivityInitial_)*(waterContentSurface_ - waterContentInitial_)
            );
        };

        return integrateScalarFunction(integrand, waterContent, waterContentRef_);
    }

    // relative transformed coordinate position (Eq. 4 Schnepf 2020)
    double deltaEtaExact(double waterContent) const
    {
        return (waterContentSurface_ - waterContentInitial_)*integral(waterContent);
    }

    // transformed coordinate
    double eta(double x /*in cm*/, double t /*in days*/) const
    {
        return std::fabs(x) - (conductivitySurface_ - conductivityInitial_)/(waterContentSurface_ - waterContentInitial_)*t;
    }

    PcKrSwCurve pcKrSwCurve_;
    double waterContentSurface_, waterContentInitial_, waterContentRef_;
    double waterContentRefDepth_, waterContentRefTime_;
    double conductivitySurface_, conductivityInitial_;
    double porosity_, permeability_, viscosity_, density_, gravity_;

    std::vector<double> waterContents_;
    std::vector<double> deltaEtas_;
};

/*!
 * \file
 * \ingroup RichardsTests
 * \brief Analytical solution for evaporation scenario
 * \note Root-soil benchmark paper Schnepf et al. (case M2.2) https://doi.org/10.3389/fpls.2020.00316
 * \note based on Vanderborght 2005 (see Fig. 5abcd and Eq. 39-47) https://doi.org/10.2113/4.1.206
 * \note Derivation in Lockington 1994 https://doi.org/10.1029/93WR03411
 *
 * The analytical solution is an approximative solution for an infinite domain
 * and is neglects graviational effects. Lockington 1994 estimates less than 1% error
 * against a given exact solution (without gravity) for a simplified pc-sw relationship.
 */
class AnalyticalSolutionM22
{
    using PcKrSwCurve = FluidMatrix::VanGenuchtenNoReg<double>;
public:
    AnalyticalSolutionM22()
    : pcKrSwCurve_("SpatialParams")
    {
        gravity_ = 9.81;
        density_ = Components::SimpleH2O<double>::liquidDensity(0,0);
        viscosity_ = Components::SimpleH2O<double>::liquidViscosity(0,0);
        porosity_ = getParam<double>("SpatialParams.Porosity");
        permeability_ = getParam<double>("SpatialParams.Permeability");

        const auto initialHead = getParam<double>("Problem.InitialHeadInCm", -40.0)*0.01; // cm to m
        potentialEvaporationRate_ = getParam<double>("Problem.SurfaceFluxMilliMeterPerDay");

        const auto pwInitial = initialHead*gravity_*density_ + 1.0e5;
        const auto pcInitial = 1.0e5 - pwInitial;
        const auto swIntial = pcKrSwCurve_.sw(pcInitial);
        waterContentInitial_ = swIntial*porosity_;

        const auto pwSurface = -10000*0.01*gravity_*density_ + 1.0e5;
        const auto pcSurface = 1.0e5 - pwSurface;
        const auto swSurface = pcKrSwCurve_.sw(pcSurface);
        waterContentSurface_ = swSurface*porosity_;

        std::cout << std::endl << "Computed analytical solution (Vanderborght 2005 5abcd, Schnepf 2020 M2.2) with:\n";
        const auto soilType = getParam<std::string>("SpatialParams.Name");
        std::cout << "Material: " << soilType << "\n";
        std::cout << Fmt::format("-> θ_i: {}\n-> θ_sur: {}\n", waterContentInitial_, waterContentSurface_);

        // water content from "effective" water content, Vanderborght 2005 Eq. 40
        const auto theta = [this](auto thetaEff){
            return thetaEff*(waterContentInitial_ - waterContentSurface_) + waterContentSurface_;
        };

        const auto dIntegral = integrateScalarFunction(
            [&,this](auto thetaEff){ return diffusivity(theta(thetaEff)); }, 0.0, 1.0
        );

        std::cout << Fmt::format("D(0), D(1): {}, {}\n", diffusivity(waterContentSurface_), diffusivity(waterContentInitial_));
        std::cout << Fmt::format("Kf(0), Kf(1): {}, {}\n", conductivity(waterContentSurface_), conductivity(waterContentInitial_));
        std::cout << Fmt::format("-> Integral D(Θ): {}\n", dIntegral);

        const auto dThetaIntegral = integrateScalarFunction(
            [&,this](auto thetaEff){ return thetaEff*diffusivity(theta(thetaEff)); }, 0.0, 1.0
        );

        std::cout << Fmt::format("0D(0), 0.5D(0.5) 1D(1): {}, {}, {}\n",
                        0.0*diffusivity(theta(0.0)), 0.5*diffusivity(theta(0.5)), 1.0*diffusivity(theta(1.0)));
        std::cout << Fmt::format("-> Integral Θ D(Θ): {}\n", dThetaIntegral);

        // Eq. 43 Vanderborght 2005
        const auto beta = dThetaIntegral*dThetaIntegral/(dIntegral*dIntegral);

        const auto dThetaFactorIntegral = integrateScalarFunction(
            [&,this](auto thetaEff){
                const auto weight = (1.0-beta*thetaEff);
                return weight*weight*diffusivity(theta(thetaEff));
            }, 0.0, 1.0
        );

        std::cout << Fmt::format("-> Integral (1.0 - βΘ)² D(θ): {}\n", dThetaFactorIntegral);

        // Eq. 42 Vanderborght 2005
        const auto alpha = dThetaFactorIntegral/dIntegral;

        // Eq. 41 Vanderborght 2005
        const auto bracketFactor = 1.0 - alpha/((1.0 - beta)*(1.0 - beta));
        const auto mu = 3*beta*(1.0 + std::sqrt(1.0 - 14.0/9.0*bracketFactor)) / ( 2*(beta - 1.0)*bracketFactor );

        std::cout << Fmt::format("-> µ: {}, α: {}, β: {}\n", mu, alpha, beta);

        // Eq. 39 Vanderborght 2005
        desorptivity_ = (waterContentInitial_ - waterContentSurface_)*std::sqrt(4.0/mu*dIntegral);

        std::cout << Fmt::format("-> Desorptivity, Sw(θ_i, θ_sur): {}\n", desorptivity_);

        // Eq. 44 & 45 Vanderborght 2005
        tPotential_ = desorptivity_*desorptivity_/(2.0*potentialEvaporationRate_*0.1*potentialEvaporationRate_*0.1);
        tPrime_ = 0.5*tPotential_;

        std::cout << Fmt::format("-> Critical time (t_pot): {} days\n", tPotential_);

        std::cout << std::endl;
    }

    // Evaporationrate in mm/day (Eq. 46 & 47 Vanderborght 2005)
    double evaporationRate(double t /*in days*/) const
    {
        if (t < tPotential_)
            return potentialEvaporationRate_;
        else
            return 0.5*desorptivity_/std::sqrt(tPrime_ + t - tPotential_)*10;
    }

private:
    // hydraulic conductivity Kf (including relative permeability) in cm/day
    double conductivity(double waterContent) const
    {
        const auto sw = waterContent/porosity_;
        const auto krw = pcKrSwCurve_.krw(sw);
        return krw*permeability_*gravity_*density_/viscosity_*100*86400;
    }

    // D = Kf*dh/dtheta Van Genuchten 1980 Eq. 10 in cm^2/day
    double diffusivity(double waterContent) const
    {
        // pc = -h*gravity_*density_
        // h = -pc/(gravity_*density_)
        const auto sw = waterContent/porosity_;
        const auto dh_dpc = -1.0/(0.01*gravity_*density_);
        const auto dpc_dsw = pcKrSwCurve_.dpc_dsw(sw);
        const auto dsw_dtheta = 1.0/porosity_;
        const auto dh_dtheta = dh_dpc*dpc_dsw*dsw_dtheta;
        return conductivity(waterContent)*dh_dtheta;
    }

    PcKrSwCurve pcKrSwCurve_;
    double waterContentSurface_, waterContentInitial_;
    double desorptivity_;
    double tPrime_, tPotential_;
    double potentialEvaporationRate_;
    double porosity_, permeability_, viscosity_, density_, gravity_;
};

} // end namespace Dumux

#endif
