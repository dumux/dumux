// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef BENCHMARKS_EXCHANGE_FLUX_CALCULATOR_HH
#define BENCHMARKS_EXCHANGE_FLUX_CALCULATOR_HH

#include <dumux/common/parameters.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>

namespace Dumux {

template <class TypeTag>
class ExchangeFluxCalculator
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    // the bulk (rock) and the low-dim (well) domain use different models, so the
    // primary variables of each domain are accessed with the indices of its own model
    using MDTraits = typename CouplingManager::MultiDomainTraits;
    using BulkTypeTag = typename MDTraits::template SubDomain<0>::TypeTag;
    using LowDimTypeTag = typename MDTraits::template SubDomain<1>::TypeTag;
    using BulkIndices = typename GetPropType<BulkTypeTag, Properties::ModelTraits>::Indices;
    using LowDimIndices = typename GetPropType<LowDimTypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<LowDimTypeTag, Properties::FluidSystem>;

    // index of the energy balance in the domain of TypeTag
    static constexpr int energyEqIdx = GetPropType<TypeTag, Properties::ModelTraits>::Indices::energyEqIdx;

    public:
    ExchangeFluxCalculator(std::shared_ptr<CouplingManager> couplingManager)
    : couplingManager_(couplingManager)
    {}

    NumEqVector computeSourceValues(const std::size_t id)
    {
        NumEqVector values(0.0);

        const Scalar temperature3D = this->couplingManager().bulkPriVars(id)[BulkIndices::temperatureIdx];
        const Scalar temperature1D = this->couplingManager().lowDimPriVars(id)[LowDimIndices::temperatureIdx];

        const Scalar radius = this->couplingManager().radius(id);
        const Scalar surface = 2*M_PI*radius;

        // Compute the convective energy flux
        values[energyEqIdx] = surface*convectionCoeff_(id)*(temperature1D - temperature3D);

        return values;
    }

    /*!
     * \brief The Nusselt number of the pipe flow.
     *
     * The flow in this benchmark is laminar (Re of about 865), so the Nusselt number of
     * laminar pipe flow with constant wall heat flux is used.
     */
    Scalar getNusseltNumber(const std::size_t id) const
    { return 4.364; }

    /*!
     * \brief The flux scaling factor of the energy exchange for the kernel coupling method.
     *
     * With the kernel method the bulk temperature entering the exchange term is an average
     * over the kernel support around the wellbore instead of the temperature at the well
     * wall. This factor rescales the exchange flux such that it matches the flux of the
     * radially symmetric analytical solution around the well, see Koch et al. (2020),
     * https://doi.org/10.1016/j.jcp.2020.109370.
     *
     * \param id The id of the point source.
     */
    Scalar energyFluxScalingFactor(const std::size_t id) const
    {
        const Scalar radius = this->couplingManager().radius(id);
        const Scalar avgDist = this->couplingManager().averageDistance(id);
        static const Scalar kernelWidthFactor = getParam<Scalar>("MixedDimension.KernelWidthFactor");
        const Scalar kernelWidth = kernelWidthFactor*radius;

        // The rock porosity is zero in this test, hence the effective thermal conductivity
        // of the bulk domain equals the thermal conductivity of the solid.
        static const Scalar lambdaBulk = getParam<Scalar>("Component.SolidThermalConductivity");

        // the dimensionless ratio of wall heat transfer coefficient to bulk conductivity
        const Scalar theta = radius*convectionCoeff_(id)/lambdaBulk;

        using std::log;
        if (avgDist <= kernelWidth)
            return 1.0/(1.0 + theta*(avgDist*avgDist/(2*kernelWidth*kernelWidth) + log(kernelWidth/radius) - 0.5));
        else
            return 1.0/(1.0 + theta*log(avgDist/radius));
    }

private:

    /*!
     * \brief Calculate the effective heat transfer coefficient including fluid convection and pipe wall conduction.
     *
     * Combines fluid-side convection (at inner surface) and cylindrical wall conduction in series,
     * referred to the outer surface (radius r_o), consistent with surface = 2*pi*r_o:
     *   1/h_eff = r_o/(r_i * h_fluid) + r_o * ln(r_o/r_i) / lambda_wall
     *
     * \param id The id of the point source.
     * \return The effective heat transfer coefficient [W/(m²·K)] referred to the outer surface.
     */
    Scalar convectionCoeff_(const std::size_t id) const
    {
        const Scalar T = this->couplingManager().lowDimPriVars(id)[LowDimIndices::temperatureIdx];
        const Scalar p = this->couplingManager().lowDimPriVars(id)[LowDimIndices::pressureIdx];
        const Scalar lambdaFluid = FluidSystem::thermalConductivity(T, p);

        const Scalar rOuter = this->couplingManager().radius(id);
        const Scalar rInner = this->couplingManager().innerRadius(id);

        // h_fluid uses inner diameter (hydraulic diameter for Nusselt correlation)
        const Scalar hFluid = lambdaFluid * getNusseltNumber(id) / (2*rInner);

        // Pipe wall conduction resistance, referred to outer surface
        static const Scalar lambdaWall = getParam<Scalar>("MixedDimension.PipeWallThermalConductivity");
        const Scalar wallResistance = rOuter * std::log(rOuter / rInner) / lambdaWall;

        // Fluid resistance mapped to outer surface: r_o/(r_i * h_fluid)
        const Scalar fluidResistance = rOuter / (rInner * hFluid);

        return 1.0 / (fluidResistance + wallResistance);
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    std::shared_ptr<CouplingManager> couplingManager_;
};


} // end namespace Dumux

#endif
