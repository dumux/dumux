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
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using MDTraits = typename CouplingManager::MultiDomainTraits;
    static constexpr auto bulkIdx = typename MDTraits::template SubDomain<0>::Index();
    static constexpr auto lowDimIdx = typename MDTraits::template SubDomain<1>::Index();

    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    public:
    ExchangeFluxCalculator(std::shared_ptr<CouplingManager> couplingManager)
    : couplingManager_(couplingManager)
    {}

    NumEqVector computeSourceValues(const std::size_t id)
    {
        NumEqVector values(0.0);

        const Scalar temperature3D = this->couplingManager().bulkPriVars(id)[temperatureIdx];
        const Scalar temperature1D = this->couplingManager().lowDimPriVars(id)[temperatureIdx];

        const Scalar radius = this->couplingManager().radius(id);
        const Scalar surface = 2*M_PI*radius;

        // Compute the convective energy flux
        values[energyEqIdx] = surface*convectionCoeff_(id)*(temperature1D - temperature3D);

        return values;
    }

    const Scalar getNusseltNumber(const std::size_t id) const
    {
        const Scalar T = this->couplingManager().lowDimPriVars(id)[temperatureIdx];
        const Scalar p = this->couplingManager().lowDimPriVars(id)[pressureIdx];
        const Scalar v = 1e-5;//this->couplingManager().lowDimMeanVelocity(id);
        const Scalar D = 2*this->couplingManager().radius(id);

        const Scalar nu =  FluidSystem::viscosity(T, p) / FluidSystem::density(T, p);
        const Scalar thermalDiffusivity = FluidSystem::thermalConductivity(T, p) / (FluidSystem::density(T, p) * FluidSystem::heatCapacity(T, p));

        const Scalar Re = abs(v) * D / nu; //Reynolds number
        const Scalar Pr = nu / thermalDiffusivity; // Prandtl number

        return NusseltNumber_(Re, Pr, D);
    }

private:

    /*!
     * \brief Calculate the convection coefficient based on the formula h = lambda_w * Nu / D
     * \param id The id of the point source.
     * \return The convection coefficient.
     */
    Scalar convectionCoeff_(const std::size_t id) const
    {
        const Scalar T = this->couplingManager().lowDimPriVars(id)[temperatureIdx];
        const Scalar p = this->couplingManager().lowDimPriVars(id)[pressureIdx];
        const Scalar lambda_w = FluidSystem::thermalConductivity(T, p);

        const Scalar D = 2*this->couplingManager().radius(id);

        return lambda_w * getNusseltNumber(id) / D;
    }

    /*!
    * \brief Calculate the Nusselt number based on the Reynolds and Prandtl numbers.
     * \param Re The Reynolds number.
     * \param Pr The Prandtl number.
     * \param D The diameter of the low-dimensional domain.
     * \return The Nusselt number.
     */
    Scalar NusseltNumber_(const Scalar Re, const Scalar Pr, const Scalar D) const
    {
        // Estimate convective film coefficient
        // const Scalar xi = 1.0 / ((1.8 * std::log(Re) - 1.5) * (1.8 * std::log(Re) - 1.5));

        // const Scalar Nu_turb = (xi / 8.0 * Re * Pr)
        //                     / (1.0 + 12.7 * std::sqrt(xi / 8.0) * (std::pow(Pr, 2.0/3.0) - 1.0))
        //                     * (1.0 + std::pow(M_PI / Z, 2.0/3.0));

        const Scalar Nu_lam = 4.364;

        return Nu_lam;

        // const Scalar gamma = (Re - 2300.0) / (10000.0 - 2300.0);

        // const Scalar Nu_m = (1.0 - gamma) * 4.364
        //                     + gamma * (0.0308 / 8.0 * 1.0e4 * Pr)
        //                     / (1.0 + 12.7 * std::sqrt(0.0308 / 8.0) * (std::pow(Pr, 2.0/3.0) - 1.0))
        //                     * (1.0 + std::pow(M_PI / Z, 2.0/3.0));

        // // Evaluate Nusselt number; based on the value of Re choose appropriate correlation
        // Scalar Nu_p;
        // if (Re < 2300.0)
        //     Nu_p = Nu_lam;
        // else if (Re > 10000.0)
        //     Nu_p = Nu_turb;
        // else
        // Nu_p = Nu_m;
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    std::shared_ptr<CouplingManager> couplingManager_;
};


} // end namespace Dumux

#endif
