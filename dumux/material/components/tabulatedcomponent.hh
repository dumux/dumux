// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Tabulates all thermodynamic properties of a given
 *        untabulated chemical species.
 *
 * At the moment, this class can only handle the sub-critical fluids
 * since it tabulates along the vapor pressure curve.
 */
#ifndef DUMUX_TABULATED_COMPONENT_HH
#define DUMUX_TABULATED_COMPONENT_HH

#include <cmath>
#include <limits>
#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <dune/common/std/type_traits.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/parallel/multithreading.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/material/components/componenttraits.hh>

namespace Dumux::Components {
// forward declaration
template<class RawComponent, bool useVaporPressure>
class TabulatedComponent;
} // end namespace Dumux::Components

namespace Dumux {

//! component traits for tabulated component
template<class RawComponent, bool useVaporPressure>
struct ComponentTraits<Components::TabulatedComponent<RawComponent, useVaporPressure>>
{
    using Scalar = typename RawComponent::Scalar;

    //! if the component implements a solid state
    static constexpr bool hasSolidState = std::is_base_of<Components::Solid<Scalar, RawComponent>, RawComponent>::value;

    //! if the component implements a liquid state
    static constexpr bool hasLiquidState = std::is_base_of<Components::Liquid<Scalar, RawComponent>, RawComponent>::value;

    //! if the component implements a gaseous state
    static constexpr bool hasGasState = std::is_base_of<Components::Gas<Scalar, RawComponent>, RawComponent>::value;
};
} // end namespace Dumux

namespace Dumux::Components::Detail {
struct DisableStaticAssert {};
} // end namespace Dumux::Components::Detail

namespace Dumux {
template<> struct AlwaysFalse<Components::Detail::DisableStaticAssert> : public std::true_type {};
}// end namespace Dumux

namespace Dumux::Components::Detail {

template<class C> using CompHasNoLiquidEnthalpy = decltype(C::template liquidEnthalpy<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoLiquidDensity = decltype(C::template liquidDensity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoLiquidThermalCond = decltype(C::template liquidThermalConductivity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoLiquidHeatCapacity = decltype(C::template liquidHeatCapacity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoLiquidViscosity = decltype(C::template liquidViscosity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasLiquidPressure = decltype(C::liquidPressure(0.0, 0.0));

template<class C> using CompHasNoGasEnthalpy = decltype(C::template gasEnthalpy<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoGasDensity = decltype(C::template gasDensity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoGasThermalCond = decltype(C::template gasThermalConductivity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoGasHeatCapacity = decltype(C::template gasHeatCapacity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasNoGasViscosity = decltype(C::template gasViscosity<DisableStaticAssert>(0.0, 0.0));
template<class C> using CompHasGasPressure = decltype(C::gasPressure(0.0, 0.0));

template<class C> constexpr inline bool hasLiquidEnthalpy()
{ return !Dune::Std::is_detected<CompHasNoLiquidEnthalpy, C>::value && ComponentTraits<C>::hasLiquidState; }
template<class C> constexpr inline bool hasLiquidDensity()
{ return !Dune::Std::is_detected<CompHasNoLiquidDensity, C>::value && ComponentTraits<C>::hasLiquidState; }
template<class C> constexpr inline bool hasLiquidThermalConductivity()
{ return !Dune::Std::is_detected<CompHasNoLiquidThermalCond, C>::value && ComponentTraits<C>::hasLiquidState; }
template<class C> constexpr inline bool hasLiquidHeatCapacity()
{ return !Dune::Std::is_detected<CompHasNoLiquidHeatCapacity, C>::value && ComponentTraits<C>::hasLiquidState; }
template<class C> constexpr inline bool hasLiquidViscosity()
{ return !Dune::Std::is_detected<CompHasNoLiquidViscosity, C>::value && ComponentTraits<C>::hasLiquidState; }
template<class C> constexpr inline bool hasLiquidPressure()
{ return Dune::Std::is_detected<CompHasLiquidPressure, C>::value && ComponentTraits<C>::hasLiquidState; }

template<class C> constexpr inline bool hasGasEnthalpy()
{ return !Dune::Std::is_detected<CompHasNoGasEnthalpy, C>::value && ComponentTraits<C>::hasGasState; }
template<class C> constexpr inline bool hasGasDensity()
{ return !Dune::Std::is_detected<CompHasNoGasDensity, C>::value && ComponentTraits<C>::hasGasState; }
template<class C> constexpr inline bool hasGasThermalConductivity()
{ return !Dune::Std::is_detected<CompHasNoGasThermalCond, C>::value && ComponentTraits<C>::hasGasState; }
template<class C> constexpr inline bool hasGasHeatCapacity()
{ return !Dune::Std::is_detected<CompHasNoGasHeatCapacity, C>::value && ComponentTraits<C>::hasGasState; }
template<class C> constexpr inline bool hasGasViscosity()
{ return !Dune::Std::is_detected<CompHasNoGasViscosity, C>::value && ComponentTraits<C>::hasGasState; }
template<class C> constexpr inline bool hasGasPressure()
{ return Dune::Std::is_detected<CompHasGasPressure, C>::value && ComponentTraits<C>::hasGasState; }

template<class RawComponent, bool useVaporPressure = true>
class TabulatedComponentTable
{
    using Scalar = typename RawComponent::Scalar;
    friend class TabulatedComponent<RawComponent, useVaporPressure>;

    struct GasPolicy
    {
        Scalar minP(std::size_t iT) const { return table.minGasPressure(iT); }
        Scalar maxP(std::size_t iT) const { return table.maxGasPressure(iT); }
        const TabulatedComponentTable<RawComponent, useVaporPressure>& table;
    };

    struct LiquidPolicy
    {
        Scalar minP(std::size_t iT) const { return table.minLiquidPressure(iT); }
        Scalar maxP(std::size_t iT) const { return table.maxLiquidPressure(iT); }
        const TabulatedComponentTable<RawComponent, useVaporPressure>& table;
    };
public:
    void init(Scalar tempMin, Scalar tempMax, std::size_t nTemp,
              Scalar pressMin, Scalar pressMax, std::size_t nPress)
    {
        tempMin_ = tempMin;
        tempMax_ = tempMax;
        nTemp_ = nTemp;
        pressMin_ = pressMin;
        pressMax_ = pressMax;
        nPress_ = nPress;
        nDensity_ = nPress;

        // resize & initialize the arrays with NaN
        assert(std::numeric_limits<Scalar>::has_quiet_NaN);
        const auto NaN = std::numeric_limits<Scalar>::quiet_NaN();

        // initialize vapor pressure array depending on useVaporPressure
        vaporPressure_.resize(nTemp_, NaN);
        tabularizeVaporPressure_();

        if constexpr (ComponentTraits<RawComponent>::hasGasState)
        {
            minGasDensity_.resize(nTemp_, NaN);
            maxGasDensity_.resize(nTemp_, NaN);
            const std::size_t numEntriesTp = nTemp_*nPress_;
            gasEnthalpy_.resize(numEntriesTp, NaN);
            gasHeatCapacity_.resize(numEntriesTp, NaN);
            gasDensity_.resize(numEntriesTp, NaN);
            gasViscosity_.resize(numEntriesTp, NaN);
            gasThermalConductivity_.resize(numEntriesTp, NaN);

            if constexpr (RawComponent::gasIsCompressible())
                gasPressure_.resize(numEntriesTp, NaN);

            minMaxGasDensityInitialized_ = tabularizeMinMaxGasDensity_();
            gasPressureInitialized_ = tabularizeGasPressure_();
            gasEnthalpyInitialized_ = tabularizeGasEnthalpy_();
            gasHeatCapacityInitialized_ = tabularizeGasHeatCapacity_();
            gasDensityInitialized_ = tabularizeGasDensity_();
            gasViscosityInitialized_ = tabularizeGasViscosity_();
            gasThermalConductivityInitialized_ = tabularizeGasThermalConductivity_();
        }

        if constexpr (ComponentTraits<RawComponent>::hasLiquidState)
        {
            minLiquidDensity_.resize(nTemp_, NaN);
            maxLiquidDensity_.resize(nTemp_, NaN);

            const std::size_t numEntriesTp = nTemp_*nPress_;
            liquidEnthalpy_.resize(numEntriesTp, NaN);
            liquidHeatCapacity_.resize(numEntriesTp, NaN);
            liquidDensity_.resize(numEntriesTp, NaN);
            liquidViscosity_.resize(numEntriesTp, NaN);
            liquidThermalConductivity_.resize(numEntriesTp, NaN);

            if constexpr (RawComponent::liquidIsCompressible())
                liquidPressure_.resize(numEntriesTp, NaN);

            minMaxLiquidDensityInitialized_ = tabularizeMinMaxLiquidDensity_();
            liquidPressureInitialized_ = tabularizeLiquidPressure_();
            liquidEnthalpyInitialized_ = tabularizeLiquidEnthalpy_();
            liquidHeatCapacityInitialized_ = tabularizeLiquidHeatCapacity_();
            liquidDensityInitialized_ = tabularizeLiquidDensity_();
            liquidViscosityInitialized_ = tabularizeLiquidViscosity_();
            liquidThermalConductivityInitialized_ = tabularizeLiquidThermalConductivity_();
        }
    }

    //! returns the index of an entry in a temperature field
    inline Scalar tempIdx(Scalar temperature) const
    {
        return (nTemp_ - 1)*(temperature - tempMin_)/(tempMax_ - tempMin_);
    }

    //! returns the index of an entry in a pressure field
    inline Scalar pressLiquidIdx(Scalar pressure, std::size_t tempIdx) const
    {
        const Scalar plMin = minLiquidPressure(tempIdx);
        const Scalar plMax = maxLiquidPressure(tempIdx);
        return (nPress_ - 1)*(pressure - plMin)/(plMax - plMin);
    }

    //! returns the index of an entry in a temperature field
    inline Scalar pressGasIdx(Scalar pressure, std::size_t tempIdx) const
    {
        const Scalar pgMin = minGasPressure(tempIdx);
        const Scalar pgMax = maxGasPressure(tempIdx);
        return (nPress_ - 1)*(pressure - pgMin)/(pgMax - pgMin);
    }

    //! returns the index of an entry in a density field
    inline Scalar densityLiquidIdx(Scalar density, std::size_t tempIdx) const
    {
        const Scalar densityMin = minLiquidDensity_[tempIdx];
        const Scalar densityMax = maxLiquidDensity_[tempIdx];
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    //! returns the index of an entry in a density field
    inline Scalar densityGasIdx(Scalar density, std::size_t tempIdx) const
    {
        const Scalar densityMin = minGasDensity_[tempIdx];
        const Scalar densityMax = maxGasDensity_[tempIdx];
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    //! returns the minimum tabularized liquid pressure at a given temperature index
    inline Scalar minLiquidPressure(int tempIdx) const
    {
        using std::max;
        if (!useVaporPressure)
            return pressMin_;
        else
            return max(pressMin_, vaporPressure_[tempIdx] / 1.1);
    }

    //! returns the maximum tabularized liquid pressure at a given temperature index
    inline Scalar maxLiquidPressure(int tempIdx) const
    {
        using std::max;
        if (!useVaporPressure)
            return pressMax_;
        else
            return max(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }

    //! returns the minimum tabularized gas pressure at a given temperature index
    inline Scalar minGasPressure(int tempIdx) const
    {
        using std::min;
        if (!useVaporPressure)
            return pressMin_;
        else
            return min(pressMin_, vaporPressure_[tempIdx] / 1.1 );
    }

    //! returns the maximum tabularized gas pressure at a given temperature index
    inline Scalar maxGasPressure(int tempIdx) const
    {
        using std::min;
        if (!useVaporPressure)
            return pressMax_;
        else
            return min(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }

    //! returns the minimum tabularized liquid density at a given temperature index
    inline Scalar minLiquidDensity(int tempIdx) const
    { return minLiquidDensity_[tempIdx]; }

    //! returns the maximum tabularized liquid density at a given temperature index
    inline Scalar maxLiquidDensity(int tempIdx) const
    { return maxLiquidDensity_[tempIdx]; }

    //! returns the minimum tabularized gas density at a given temperature index
    inline Scalar minGasDensity(int tempIdx) const
    { return minGasDensity_[tempIdx]; }

    //! returns the maximum tabularized gas density at a given temperature index
    inline Scalar maxGasDensity(int tempIdx) const
    { return maxGasDensity_[tempIdx]; }

    inline std::size_t nTemp() const { return nTemp_; }
    inline std::size_t nPress() const { return nPress_; }
    inline std::size_t nDensity() const { return nDensity_; }

    inline Scalar tempMax() const { return tempMax_; }
    inline Scalar tempMin() const { return tempMin_; }

private:

    //! initializes vapor pressure if useVaporPressure = true
    template< bool useVP = useVaporPressure, std::enable_if_t<useVP, int> = 0 >
    void tabularizeVaporPressure_()
    {
        // fill the temperature-pressure arrays
        Dumux::parallelFor(nTemp_, [=](std::size_t iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;
            vaporPressure_[iT] = RawComponent::vaporPressure(temperature);
        });
    }

    //! if !useVaporPressure, do nothing here
    template< bool useVP = useVaporPressure, std::enable_if_t<!useVP, int> = 0 >
    void tabularizeVaporPressure_() {}

    /*!
     * \brief Initializes property values as function of temperature and pressure.
     *
     * \tparam PropFunc Function to evaluate the property prop(T, p)
     * \tparam MinPFunc Function to evaluate the minimum pressure for a
     *                  temperature index (depends on useVaporPressure)
     * \tparam MaxPFunc Function to evaluate the maximum pressure for a
     *                  temperature index (depends on useVaporPressure)
     *
     * \param f property function
     * \param minP function to evaluate minimum pressure for temp idx
     * \param maxP function to evaluate maximum pressure for temp idx
     * \param values container to store property values
     */
    template<class PropFunc, class Policy>
    void tabularizeTPArray_(const PropFunc& f, Policy policy, std::vector<typename RawComponent::Scalar>& values) const
    {
        Dumux::parallelFor(nTemp_, [=,&values](std::size_t iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            Scalar pMax = policy.maxP(iT);
            Scalar pMin = policy.minP(iT);
            for (std::size_t iP = 0; iP < nPress_; ++ iP)
            {
                Scalar pressure = iP * (pMax - pMin)/(nPress_ - 1) + pMin;
                values[iT + iP*nTemp_] = f(temperature, pressure);
            }
        });
    }

    /*!
     * \brief Initializes the minimum/maximum densities on the temperature range.
     *
     * \tparam RhoFunc Function to evaluate the density rho(T, p)
     * \tparam MinPFunc Function to evaluate the minimum pressure for a
     *                  temperature index (depends on useVaporPressure)
     * \tparam MaxPFunc Function to evaluate the maximum pressure for a
     *                  temperature index (depends on useVaporPressure)
     *
     * \param rho density function
     * \param minP function to evaluate minimum pressure for temp idx
     * \param maxP function to evaluate maximum pressure for temp idx
     * \param rhoMin container to store minimum density values
     * \param rhoMax container to store maximum density values
     */
    template<class RhoFunc, class Policy>
    void tabularizeMinMaxRhoArray_(const RhoFunc& rho, Policy policy,
                                   std::vector<typename RawComponent::Scalar>& rhoMin,
                                   std::vector<typename RawComponent::Scalar>& rhoMax) const
    {
        Dumux::parallelFor(nTemp_, [=,&rhoMin,&rhoMax](std::size_t iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            using std::min;
            rhoMin[iT] = rho(temperature, policy.minP(iT));
            rhoMax[iT] = rho(temperature, policy.maxP(min(iT + 1, nTemp_ - 1)));
        });
    }

    /*!
     * \brief Initializes pressure arrays as function of temperature and density.
     *
     * \tparam PFunc Function to evaluate the pressure p(T, rho)
     *
     * \param pressure container to store pressure values
     * \param p pressure function p(T, rho)
     * \param rhoMin container with minimum density values
     * \param rhoMax container with maximum density values
     */
    template<class PFunc>
    void tabularizePressureArray_(std::vector<typename RawComponent::Scalar>& pressure,
                                  const PFunc& p,
                                  const std::vector<typename RawComponent::Scalar>& rhoMin,
                                  const std::vector<typename RawComponent::Scalar>& rhoMax) const
    {
        Dumux::parallelFor(nTemp_, [=,&pressure,&rhoMin,&rhoMax](std::size_t iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            for (std::size_t iRho = 0; iRho < nDensity_; ++ iRho)
            {
                Scalar density = Scalar(iRho)/(nDensity_ - 1)
                                 * (rhoMax[iT] - rhoMin[iT])
                                 +  rhoMin[iT];
                pressure[iT + iRho*nTemp_] = p(temperature, density);
            }
        });
    }

    template<class RC = RawComponent>
    bool tabularizeGasEnthalpy_()
    {
        if constexpr (Detail::hasGasEnthalpy<RC>())
        {
            const auto gasEnth = [] (auto T, auto p) { return RC::gasEnthalpy(T, p); };
            tabularizeTPArray_(gasEnth, GasPolicy{ *this }, gasEnthalpy_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeLiquidEnthalpy_()
    {
        if constexpr (Detail::hasLiquidEnthalpy<RC>())
        {
            const auto liqEnth = [] (auto T, auto p) { return RC::liquidEnthalpy(T, p); };
            tabularizeTPArray_(liqEnth, LiquidPolicy{ *this }, liquidEnthalpy_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeGasHeatCapacity_()
    {
        if constexpr (Detail::hasGasHeatCapacity<RC>())
        {
            const auto gasHC = [] (auto T, auto p) { return RC::gasHeatCapacity(T, p); };
            tabularizeTPArray_(gasHC, GasPolicy{ *this }, gasHeatCapacity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeLiquidHeatCapacity_()
    {
        if constexpr (Detail::hasLiquidHeatCapacity<RC>())
        {
            const auto liqHC = [] (auto T, auto p) { return RC::liquidHeatCapacity(T, p); };
            tabularizeTPArray_(liqHC, LiquidPolicy{ *this }, liquidHeatCapacity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeMinMaxGasDensity_()
    {
        if constexpr (Detail::hasGasDensity<RC>())
        {
            const auto gasRho = [] (auto T, auto p) { return RC::gasDensity(T, p); };
            tabularizeMinMaxRhoArray_(gasRho, GasPolicy{ *this }, minGasDensity_, maxGasDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeMinMaxLiquidDensity_()
    {
        if constexpr (Detail::hasGasEnthalpy<RC>())
        {
            const auto liqRho = [] (auto T, auto p) { return RC::liquidDensity(T, p); };
            tabularizeMinMaxRhoArray_(liqRho, LiquidPolicy{ *this }, minLiquidDensity_, maxLiquidDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeGasPressure_()
    {
        // pressure is only defined if the gas is compressible (this is usually the case)
        if constexpr (Detail::hasGasPressure<RC>() && RC::gasIsCompressible())
        {
            const auto gasPFunc = [] (auto T, auto rho) { return RC::gasPressure(T, rho); };
            tabularizePressureArray_(gasPressure_, gasPFunc, minGasDensity_, maxGasDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeLiquidPressure_()
    {
        // pressure is only defined if the liquid is compressible (this is often not the case)
        if constexpr (Detail::hasLiquidPressure<RC>() && RC::liquidIsCompressible())
        {
            const auto liqPFunc = [] (auto T, auto rho) { return RC::liquidPressure(T, rho); };
            tabularizePressureArray_(liquidPressure_, liqPFunc, minLiquidDensity_, maxLiquidDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeGasDensity_()
    {
        if constexpr (Detail::hasGasDensity<RC>())
        {
            const auto gasRho = [] (auto T, auto p) { return RC::gasDensity(T, p); };
            tabularizeTPArray_(gasRho, GasPolicy{ *this }, gasDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeLiquidDensity_()
    {
        if constexpr (Detail::hasLiquidDensity<RC>())
        {
            // TODO: we could get rid of the lambdas and pass the functor directly. But,
            //       currently Brine is a component (and not a fluid system) expecting a
            //       third argument with a default, which cannot be wrapped in a function pointer.
            //       For this reason we have to wrap this into a lambda here.
            const auto liqRho = [] (auto T, auto p) { return RC::liquidDensity(T, p); };
            tabularizeTPArray_(liqRho, LiquidPolicy{ *this }, liquidDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeGasViscosity_()
    {
        if constexpr (Detail::hasGasViscosity<RC>())
        {
            const auto gasVisc = [] (auto T, auto p) { return RC::gasViscosity(T, p); };
            tabularizeTPArray_(gasVisc, GasPolicy{ *this }, gasViscosity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeLiquidViscosity_()
    {
        if constexpr (Detail::hasLiquidViscosity<RC>())
        {
            const auto liqVisc = [] (auto T, auto p) { return RC::liquidViscosity(T, p); };
            tabularizeTPArray_(liqVisc, LiquidPolicy{ *this }, liquidViscosity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeGasThermalConductivity_()
    {
        if constexpr (Detail::hasGasThermalConductivity<RC>())
        {
            const auto gasTC = [] (auto T, auto p) { return RC::gasThermalConductivity(T, p); };
            tabularizeTPArray_(gasTC, GasPolicy{ *this }, gasThermalConductivity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    bool tabularizeLiquidThermalConductivity_()
    {
        if constexpr (Detail::hasLiquidThermalConductivity<RC>())
        {
            const auto liqTC = [] (auto T, auto p) { return RC::liquidThermalConductivity(T, p); };
            tabularizeTPArray_(liqTC, LiquidPolicy{ *this }, liquidThermalConductivity_);
            return true;
        }

        return false;
    }

private:
    // 1D fields with the temperature as degree of freedom
    std::vector<Scalar> vaporPressure_;

    std::vector<Scalar> minLiquidDensity_;
    std::vector<Scalar> maxLiquidDensity_;
    bool minMaxLiquidDensityInitialized_;

    std::vector<Scalar> minGasDensity_;
    std::vector<Scalar> maxGasDensity_;
    bool minMaxGasDensityInitialized_;

    // 2D fields with the temperature and pressure as degrees of freedom
    std::vector<Scalar> gasEnthalpy_;
    std::vector<Scalar> liquidEnthalpy_;
    bool gasEnthalpyInitialized_;
    bool liquidEnthalpyInitialized_;

    std::vector<Scalar> gasHeatCapacity_;
    std::vector<Scalar> liquidHeatCapacity_;
    bool gasHeatCapacityInitialized_;
    bool liquidHeatCapacityInitialized_;

    std::vector<Scalar> gasDensity_;
    std::vector<Scalar> liquidDensity_;
    bool gasDensityInitialized_;
    bool liquidDensityInitialized_;

    std::vector<Scalar> gasViscosity_;
    std::vector<Scalar> liquidViscosity_;
    bool gasViscosityInitialized_;
    bool liquidViscosityInitialized_;

    std::vector<Scalar> gasThermalConductivity_;
    std::vector<Scalar> liquidThermalConductivity_;
    bool gasThermalConductivityInitialized_;
    bool liquidThermalConductivityInitialized_;

    // 2D fields with the temperature and density as degrees of freedom
    std::vector<Scalar> gasPressure_;
    std::vector<Scalar> liquidPressure_;
    bool gasPressureInitialized_;
    bool liquidPressureInitialized_;

    // temperature, pressure and density ranges
    Scalar tempMin_;
    Scalar tempMax_;
    std::size_t nTemp_;

    Scalar pressMin_;
    Scalar pressMax_;
    std::size_t nPress_;

    Scalar densityMin_;
    Scalar densityMax_;
    std::size_t nDensity_;
};

} // end namespace Dumux::Components::Detail

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief  Tabulates all thermodynamic properties of a given component
 *
 * At the moment, this class can only handle the sub-critical fluids
 * since it tabulates along the vapor pressure curve.
 *
 * \tparam Scalar The type used for scalar values
 * \tparam RawComponent The component which ought to be tabulated
 * \tparam useVaporPressure If set to true, the min/max pressure
 *                          values for gas&liquid phase will be set
 *                          depending on the vapor pressure.
 */
template <class RawComponent, bool useVaporPressure=true>
class TabulatedComponent
{
    using ThisType = TabulatedComponent<RawComponent, useVaporPressure>;
    using Table = Detail::TabulatedComponentTable<RawComponent, useVaporPressure>;

    struct InterpolatePolicy
    {
        using Scalar = typename RawComponent::Scalar;

        const Table& table;

        Scalar tempIdx(Scalar T) const { return table.tempIdx(T); }
        Scalar tempMin() const { return table.tempMin(); }
        Scalar tempMax() const { return table.tempMax(); }
        std::size_t nTemp() const { return table.nTemp(); }
        std::size_t nPress() const { return table.nPress(); }
        std::size_t nDensity() const { return table.nDensity(); }
    };

    struct InterpolateGasPolicy : public InterpolatePolicy
    {
        using Scalar = typename RawComponent::Scalar;
        Scalar pressIdx(Scalar p, std::size_t tempIdx) const { return this->table.pressGasIdx(p, tempIdx); }
        Scalar rhoIdx(Scalar rho, std::size_t tempIdx) const { return this->table.densityGasIdx(rho, tempIdx); }
        Scalar minP(int tempIdx) const { return this->table.minGasPressure(tempIdx); }
        Scalar maxP(int tempIdx) const { return this->table.maxGasPressure(tempIdx); }
        Scalar minRho(int tempIdx) const { return this->table.minGasDensity(tempIdx); }
        Scalar maxRho(int tempIdx) const { return this->table.maxGasDensity(tempIdx); }
    };

    struct InterpolateLiquidPolicy : public InterpolatePolicy
    {
        using Scalar = typename RawComponent::Scalar;
        Scalar pressIdx(Scalar p, std::size_t tempIdx) const { return this->table.pressLiquidIdx(p, tempIdx); }
        Scalar rhoIdx(Scalar rho, std::size_t tempIdx) const { return this->table.densityLiquidIdx(rho, tempIdx); }
        Scalar minP(int tempIdx) const { return this->table.minLiquidPressure(tempIdx); }
        Scalar maxP(int tempIdx) const { return this->table.maxLiquidPressure(tempIdx); }
        Scalar minRho(int tempIdx) const { return this->table.minLiquidDensity(tempIdx); }
        Scalar maxRho(int tempIdx) const { return this->table.maxLiquidDensity(tempIdx); }
    };
public:
    //! export scalar type
    using Scalar = typename RawComponent::Scalar;

    //! state that we are tabulated
    static constexpr bool isTabulated = true;

    /*!
     * \brief Initialize the tables.
     *
     * \param tempMin The minimum of the temperature range in \f$\mathrm{[K]}\f$
     * \param tempMax The maximum of the temperature range in \f$\mathrm{[K]}\f$
     * \param nTemp The number of entries/steps within the temperature range
     * \param pressMin The minimum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param nPress The number of entries/steps within the pressure range
     */
    static void init(Scalar tempMin, Scalar tempMax, std::size_t nTemp,
                     Scalar pressMin, Scalar pressMax, std::size_t nPress)
    {
#ifndef NDEBUG
        warningPrinted_ = false;
#endif
        std::cout << "-------------------------------------------------------------------------\n"
                  << "Initializing tables for the " << RawComponent::name()
                  << " fluid properties (" << nTemp*nPress  << " entries).\n"
                  << "Temperature -> min: " << std::scientific << std::setprecision(3)
                  << tempMin << ", max: "  << tempMax << ", n: " << nTemp << '\n'
                  << "Pressure    -> min: " << std::scientific << std::setprecision(3)
                  << pressMin << ", max: " << pressMax << ", n: " << nPress << '\n'
                  << "-------------------------------------------------------------------------" << std::endl;

        table().init(tempMin, tempMax, nTemp, pressMin, pressMax, nPress);

        initialized_  = true;
    }

    /*!
     * \brief A human readable name for the component.
     */
    static std::string name()
    { return RawComponent::name(); }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    static constexpr Scalar molarMass()
    { return RawComponent::molarMass(); }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component.
     */
    static Scalar criticalTemperature()
    { return RawComponent::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component.
     */
    static Scalar criticalPressure()
    { return RawComponent::criticalPressure(); }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature()
    { return RawComponent::tripleTemperature(); }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure()
    { return RawComponent::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     *
     * \param T temperature of component
     */
    static Scalar vaporPressure(Scalar T)
    {
        using std::isnan;
        Scalar result = interpolateT_(table().vaporPressure_, T, InterpolatePolicy{{ table() }});
        if (isnan(result))
            return RawComponent::vaporPressure(T);
        return result;
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     *
     * The method is only called by the sequential flash, so tabulating is omitted.
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar vaporTemperature(Scalar pressure)
    {
        return RawComponent::vaporTemperature(pressure);
    }

    /*!
     * \brief Specific enthalpy of the gas \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().gasEnthalpy_, temperature, pressure, InterpolateGasPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("gasEnthalpy", temperature, pressure, InterpolateGasPolicy{{ table() }});
            return RawComponent::gasEnthalpy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific enthalpy of the liquid \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().liquidEnthalpy_, temperature, pressure, InterpolateLiquidPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("liquidEnthalpy", temperature, pressure, InterpolateLiquidPolicy{{ table() }});
            return RawComponent::liquidEnthalpy(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific isobaric heat capacity of the gas \f$\mathrm{[J/(kg*K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().gasHeatCapacity_, temperature, pressure, InterpolateGasPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("gasHeatCapacity", temperature, pressure, InterpolateGasPolicy{{ table() }});
            return RawComponent::gasHeatCapacity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific isobaric heat capacity of the liquid \f$\mathrm{[J/(kg*K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().liquidHeatCapacity_, temperature, pressure, InterpolateLiquidPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("liquidHeatCapacity", temperature, pressure, InterpolateLiquidPolicy{{ table() }});
            return RawComponent::liquidHeatCapacity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief Specific internal energy of the gas \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        return gasEnthalpy(temperature, pressure) - pressure/gasDensity(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of the liquid \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    {
        return liquidEnthalpy(temperature, pressure) - pressure/liquidDensity(temperature, pressure);
    }

    /*!
     * \brief The pressure of gas in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        Scalar result = interpolateTRho_(table().gasPressure_, temperature, density, InterpolateGasPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTRho_("gasPressure", temperature, density, InterpolateGasPolicy{{ table() }});
            return RawComponent::gasPressure(temperature, density);
        }
        return result;
    }

    /*!
     * \brief The pressure of liquid in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        Scalar result = interpolateTRho_(table().liquidPressure_, temperature, density, InterpolateLiquidPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTRho_("liquidPressure", temperature, density, InterpolateLiquidPolicy{{ table() }});
            return RawComponent::liquidPressure(temperature, density);
        }
        return result;
    }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return RawComponent::gasIsCompressible(); }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return RawComponent::liquidIsCompressible(); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return RawComponent::gasIsIdeal(); }


    /*!
     * \brief The density of gas at a given pressure and temperature
     *        \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().gasDensity_, temperature, pressure, InterpolateGasPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("gasDensity", temperature, pressure, InterpolateGasPolicy{{ table() }});
            return RawComponent::gasDensity(temperature, pressure);
        }
        return result;
    }

    /*!
     *  \brief The molar density of gas in \f$\mathrm{[mol/m^3]}\f$
     *  at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return gasDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The density of liquid at a given pressure and
     *        temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().liquidDensity_, temperature, pressure, InterpolateLiquidPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("liquidDensity", temperature, pressure, InterpolateLiquidPolicy{{ table() }});
            return RawComponent::liquidDensity(temperature, pressure);
        }

        return result;
    }

    /*!
     * \brief The molar density of liquid in \f$\mathrm{[mol/m^3]}\f$
     * at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().gasViscosity_, temperature, pressure, InterpolateGasPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("gasViscosity", temperature, pressure, InterpolateGasPolicy{{ table() }});
            return RawComponent::gasViscosity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().liquidViscosity_, temperature, pressure, InterpolateLiquidPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("liquidViscosity",temperature, pressure, InterpolateLiquidPolicy{{ table() }});
            return RawComponent::liquidViscosity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The thermal conductivity of gaseous water \f$\mathrm{[W/(m*K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().gasThermalConductivity_, temperature, pressure, InterpolateGasPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("gasThermalConductivity", temperature, pressure, InterpolateGasPolicy{{ table() }});
            return RawComponent::gasThermalConductivity(temperature, pressure);
        }
        return result;
    }

    /*!
     * \brief The thermal conductivity of liquid water \f$\mathrm{[W/(m*K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    {
        Scalar result = interpolateTP_(table().liquidThermalConductivity_, temperature, pressure, InterpolateLiquidPolicy{{ table() }});
        using std::isnan;
        if (isnan(result))
        {
            printWarningTP_("liquidThermalConductivity", temperature, pressure, InterpolateLiquidPolicy{{ table() }});
            return RawComponent::liquidThermalConductivity(temperature, pressure);
        }
        return result;
    }


private:
    //! prints a warning if the result is not in range or the table has not been initialized
    template<class Policy>
    static void printWarningTP_(const std::string& quantity, Scalar T, Scalar p, Policy policy)
    {
#ifndef NDEBUG
        if (warningPrinted_)
            return;

        if (!initialized_)
            std::cerr << "Warning: tabulated component '" << name()
                      << "' has not been initialized. "
                      << "Call FluidSystem::init() to use the tabulation in order to reduce runtime. \n";
        else
            std::cerr << "Warning: "<<quantity<<"(T="<<T<<", p="<<p<<") of component '"<<name()
                      << "' is outside tabulation range: ("<< policy.tempMin() <<"<=T<="<< policy.tempMax() <<"), ("
                      << policy.minP(0) <<"<=p<=" << policy.maxP(policy.nTemp()-1) <<"). "
                      << "Forwarded to FluidSystem for direct evaluation of "<<quantity<<". \n";
        warningPrinted_ = true;
#endif
    }

    //! prints a warning if the result is not in range or the table has not been initialized
    template<class Policy>
    static void printWarningTRho_(const std::string& quantity, Scalar T, Scalar rho, Policy policy)
    {
#ifndef NDEBUG
        if (warningPrinted_)
            return;

        if (!initialized_)
            std::cerr << "Warning: tabulated component '" << name()
                      << "' has not been initialized. "
                      << "Call FluidSystem::init() to use the tabulation in order to reduce runtime. \n";
        else
        {
            const auto [densityMin, densityMax] = [&]
            {
                const Scalar alphaT = policy.tempIdx(T);
                using std::clamp;
                const auto iT = clamp<int>(static_cast<int>(alphaT), 0, policy.nTemp() - 2);
                return std::make_pair( policy.minRho(iT), policy.maxRho(iT) );
            }();

            std::cerr << "Warning: "<<quantity<<"(T="<<T<<", density="<<rho<<") of component '"<<name()
                      << "' is outside tabulation range: ("<< policy.tempMin() <<"<=T<="<< policy.tempMax() <<"), ("
                      << densityMin <<"<=density<=" << densityMax <<"). "
                      << "Forwarded to FluidSystem for direct evaluation of "<<quantity<<". \n";
        }
        warningPrinted_ = true;
#endif
    }

    //! returns an interpolated value depending on temperature
    template<class Policy>
    static Scalar interpolateT_(const std::vector<Scalar>& values, Scalar T, Policy policy)
    {
        Scalar alphaT = policy.tempIdx(T);
        const auto nTemp = policy.nTemp();
        if (alphaT < 0 - 1e-7*nTemp || alphaT >= nTemp - 1 + 1e-7*nTemp)
            return std::numeric_limits<Scalar>::quiet_NaN();

        using std::clamp;
        const auto iT = clamp<int>(static_cast<int>(alphaT), 0, nTemp - 2);
        alphaT -= iT;

        return values[iT    ]*(1 - alphaT) +
               values[iT + 1]*(    alphaT);
    }

    //! returns an interpolated value depending on temperature and pressure
    template<class Policy>
    static Scalar interpolateTP_(const std::vector<Scalar>& values, const Scalar T, const Scalar p, Policy policy)
    {
        const auto nTemp = policy.nTemp();
        const auto nPress = policy.nPress();

        Scalar alphaT = policy.tempIdx(T);
        if (alphaT < 0 - 1e-7*nTemp || alphaT >= nTemp - 1 + 1e-7*nTemp)
            return std::numeric_limits<Scalar>::quiet_NaN();

        using std::clamp;
        const auto iT = clamp<int>(static_cast<int>(alphaT), 0, nTemp - 2);
        alphaT -= iT;

        Scalar alphaP1 = policy.pressIdx(p, iT);
        Scalar alphaP2 = policy.pressIdx(p, iT + 1);

        const auto iP1 = clamp<int>(static_cast<int>(alphaP1), 0, nPress - 2);
        const auto iP2 = clamp<int>(static_cast<int>(alphaP2), 0, nPress - 2);
        alphaP1 -= iP1;
        alphaP2 -= iP2;

#if 0 && !defined NDEBUG
        const auto tempMin = policy.tempMin();
        const auto tempMin = policy.tempMax();
        if(!(0 <= alphaT && alphaT <= 1.0))
            DUNE_THROW(NumericalProblem, "Temperature out of range: "
                       << "T=" << T << " range: [" << tempMin_ << ", " << tempMax_ << "]");
        if(!(0 <= alphaP1 && alphaP1 <= 1.0))
            DUNE_THROW(NumericalProblem, "First pressure out of range: "
                       << "p=" << p << " range: [" << policy.minP(policy.tempIdx(T)) << ", " << policy.maxP(policy.tempIdx(T)) << "]");
        if(!(0 <= alphaP2 && alphaP2 <= 1.0))
            DUNE_THROW(NumericalProblem, "Second pressure out of range: "
                       << "p=" << p << " range: [" << policy.minP(policy.tempIdx(T) + 1) << ", " << policy.maxP(policy.tempIdx(T) + 1) << "]");
#endif
        return values[(iT    ) + (iP1    )*nTemp]*(1 - alphaT)*(1 - alphaP1) +
               values[(iT    ) + (iP1 + 1)*nTemp]*(1 - alphaT)*(    alphaP1) +
               values[(iT + 1) + (iP2    )*nTemp]*(    alphaT)*(1 - alphaP2) +
               values[(iT + 1) + (iP2 + 1)*nTemp]*(    alphaT)*(    alphaP2);
    }

    //! returns an interpolated value for gas depending on temperature and density
    template<class Policy>
    static Scalar interpolateTRho_(const std::vector<Scalar>& values, const Scalar T, const Scalar rho, Policy policy)
    {
        const auto nTemp = policy.nTemp();
        const auto nDensity = policy.nDensity();

        using std::clamp;
        Scalar alphaT = policy.tempIdx(T);
        if (alphaT < 0 - 1e-7*nTemp || alphaT >= nTemp - 1 + 1e-7*nTemp)
            return std::numeric_limits<Scalar>::quiet_NaN();

        const auto iT = clamp<int>(static_cast<int>(alphaT), 0, nTemp - 2);
        alphaT -= iT;

        Scalar alphaP1 = policy.rhoIdx(rho, iT);
        Scalar alphaP2 = policy.rhoIdx(rho, iT + 1);
        const auto iP1 = clamp<int>(static_cast<int>(alphaP1), 0, nDensity - 2);
        const auto iP2 = clamp<int>(static_cast<int>(alphaP2), 0, nDensity - 2);
        alphaP1 -= iP1;
        alphaP2 -= iP2;

        return values[(iT    ) + (iP1    )*nTemp]*(1 - alphaT)*(1 - alphaP1) +
               values[(iT    ) + (iP1 + 1)*nTemp]*(1 - alphaT)*(    alphaP1) +
               values[(iT + 1) + (iP2    )*nTemp]*(    alphaT)*(1 - alphaP2) +
               values[(iT + 1) + (iP2 + 1)*nTemp]*(    alphaT)*(    alphaP2);
    }

    // specifies whether the table was initialized
    static bool initialized_;

#ifndef NDEBUG
    // specifies whether some warning was printed
    static bool warningPrinted_;
#endif

    static Table& table()
    {
        static Table t;
        return t;
    }
};

template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::initialized_ = false;

#ifndef NDEBUG
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::warningPrinted_ = false;
#endif

// forward declaration
template <class Component>
struct IsAqueous;

// we are aqueous if the raw compont is so
template <class RawComponent, bool useVaporPressure>
struct IsAqueous<TabulatedComponent<RawComponent, useVaporPressure>> : public IsAqueous<RawComponent> {};

} // end namespace Dumux::Components

#endif
