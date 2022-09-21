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
        tempMin_ = tempMin;
        tempMax_ = tempMax;
        nTemp_ = nTemp;
        pressMin_ = pressMin;
        pressMax_ = pressMax;
        nPress_ = nPress;
        nDensity_ = nPress;

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

        // resize & initialize the arrays with NaN
        assert(std::numeric_limits<Scalar>::has_quiet_NaN);
        const auto NaN = std::numeric_limits<Scalar>::quiet_NaN();

        // initialize vapor pressure array depending on useVaporPressure
        vaporPressure_.resize(nTemp_, NaN);
        initVaporPressure_();

        if constexpr (ComponentTraits<ThisType>::hasGasState)
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

            minMaxGasDensityInitialized_ = false;
            gasEnthalpyInitialized_ = false;
            gasHeatCapacityInitialized_ = false;
            gasDensityInitialized_ = false;
            gasViscosityInitialized_ = false;
            gasThermalConductivityInitialized_ = false;
            gasPressureInitialized_ = false;
        }

        if constexpr (ComponentTraits<ThisType>::hasLiquidState)
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

            // reset all flags
            minMaxLiquidDensityInitialized_ = false;
            liquidEnthalpyInitialized_ = false;
            liquidHeatCapacityInitialized_ = false;
            liquidDensityInitialized_ = false;
            liquidViscosityInitialized_ = false;
            liquidThermalConductivityInitialized_ = false;
            liquidPressureInitialized_ = false;
        }

        // in the multithreaded case precompute all array because lazy-loading
        // leads to potential data races
        if constexpr (!Multithreading::isSerial())
        {
            if constexpr (ComponentTraits<ThisType>::hasGasState)
            {
                minMaxGasDensityInitialized_ = tabularizeMinMaxGasDensity_();
                gasPressureInitialized_ = tabularizeGasPressure_();
                gasEnthalpyInitialized_ = tabularizeGasEnthalpy_();
                gasHeatCapacityInitialized_ = tabularizeGasHeatCapacity_();
                gasDensityInitialized_ = tabularizeGasDensity_();
                gasViscosityInitialized_ = tabularizeGasViscosity_();
                gasThermalConductivityInitialized_ = tabularizeGasThermalConductivity_();
            }

            if constexpr (ComponentTraits<ThisType>::hasLiquidState)
            {
                minMaxLiquidDensityInitialized_ = tabularizeMinMaxLiquidDensity_();
                liquidPressureInitialized_ = tabularizeLiquidPressure_();
                liquidEnthalpyInitialized_ = tabularizeLiquidEnthalpy_();
                liquidHeatCapacityInitialized_ = tabularizeLiquidHeatCapacity_();
                liquidDensityInitialized_ = tabularizeLiquidDensity_();
                liquidViscosityInitialized_ = tabularizeLiquidViscosity_();
                liquidThermalConductivityInitialized_ = tabularizeLiquidThermalConductivity_();
            }
        }

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
        Scalar result = interpolateT_(vaporPressure_, T);
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
        Scalar result = interpolateTP_(gasEnthalpy_, temperature, pressure,
                                       pressGasIdx_, minGasPressure_, maxGasPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!gasEnthalpyInitialized_)
                {
                    gasEnthalpyInitialized_ = tabularizeGasEnthalpy_();
                    if (gasEnthalpyInitialized_)
                        return gasEnthalpy(temperature, pressure);
                }
            }

            printWarning_("gasEnthalpy", temperature, pressure);
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
        Scalar result = interpolateTP_(liquidEnthalpy_, temperature, pressure,
                                       pressLiquidIdx_, minLiquidPressure_, maxLiquidPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!liquidEnthalpyInitialized_)
                {
                    liquidEnthalpyInitialized_ = tabularizeLiquidEnthalpy_();
                    if (liquidEnthalpyInitialized_)
                        return liquidEnthalpy(temperature, pressure);
                }
            }

            printWarning_("liquidEnthalpy", temperature, pressure);
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
        Scalar result = interpolateTP_(gasHeatCapacity_, temperature, pressure,
                                       pressGasIdx_, minGasPressure_, maxGasPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!gasHeatCapacityInitialized_)
                {
                    gasHeatCapacityInitialized_ = tabularizeGasHeatCapacity_();
                    if (gasHeatCapacityInitialized_)
                        return gasHeatCapacity(temperature, pressure);
                }
            }

            printWarning_("gasHeatCapacity", temperature, pressure);
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
        Scalar result = interpolateTP_(liquidHeatCapacity_, temperature, pressure,
                                       pressLiquidIdx_, minLiquidPressure_, maxLiquidPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!liquidHeatCapacityInitialized_)
                {
                    liquidHeatCapacityInitialized_ = tabularizeLiquidHeatCapacity_();
                    if (liquidHeatCapacityInitialized_)
                        return liquidHeatCapacity(temperature, pressure);
                }
            }

            printWarning_("liquidHeatCapacity", temperature, pressure);
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
        if constexpr (Multithreading::isSerial())
        {
            if (!minMaxGasDensityInitialized_)
                minMaxGasDensityInitialized_ = tabularizeMinMaxGasDensity_();
        }

        Scalar result = interpolateTRho_(gasPressure_, temperature, density, densityGasIdx_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                if (!gasPressureInitialized_)
                {
                    gasPressureInitialized_ = tabularizeGasPressure_();
                    if (gasPressureInitialized_)
                        return gasPressure(temperature, density);
                }
            }

            printWarning_("gasPressure", temperature, density);
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
        if constexpr (Multithreading::isSerial())
        {
            if (!minMaxLiquidDensityInitialized_)
                minMaxLiquidDensityInitialized_ = tabularizeMinMaxLiquidDensity_();
        }

        Scalar result = interpolateTRho_(liquidPressure_, temperature, density, densityLiquidIdx_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!liquidPressureInitialized_)
                {
                    liquidPressureInitialized_ = tabularizeLiquidPressure_();
                    if (liquidPressureInitialized_)
                        return liquidPressure(temperature, density);
                }
            }

            printWarning_("liquidPressure", temperature, density);
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
        Scalar result = interpolateTP_(gasDensity_, temperature, pressure,
                                       pressGasIdx_, minGasPressure_, maxGasPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!gasDensityInitialized_)
                {
                    gasDensityInitialized_ = tabularizeGasDensity_();
                    if (gasDensityInitialized_)
                        return gasDensity(temperature, pressure);
                }
            }

            printWarning_("gasDensity", temperature, pressure);
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
        Scalar result = interpolateTP_(liquidDensity_, temperature, pressure,
                                       pressLiquidIdx_, minLiquidPressure_, maxLiquidPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!liquidDensityInitialized_)
                {
                    liquidDensityInitialized_ = tabularizeLiquidDensity_();
                    if (liquidDensityInitialized_)
                        return liquidDensity(temperature, pressure);
                }
            }

            printWarning_("liquidDensity", temperature, pressure);
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
        Scalar result = interpolateTP_(gasViscosity_, temperature, pressure,
                                       pressGasIdx_, minGasPressure_, maxGasPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!gasViscosityInitialized_)
                {
                    gasViscosityInitialized_ = tabularizeGasViscosity_();
                    if (gasViscosityInitialized_)
                        return gasViscosity(temperature, pressure);
                }
            }

            printWarning_("gasViscosity", temperature, pressure);
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
        Scalar result = interpolateTP_(liquidViscosity_, temperature, pressure,
                                       pressLiquidIdx_, minLiquidPressure_, maxLiquidPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!liquidViscosityInitialized_)
                {
                    liquidViscosityInitialized_ = tabularizeLiquidViscosity_();
                    if (liquidViscosityInitialized_)
                        return liquidViscosity(temperature, pressure);
                }
            }

            printWarning_("liquidViscosity",temperature, pressure);
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
        Scalar result = interpolateTP_(gasThermalConductivity_, temperature, pressure,
                                       pressGasIdx_, minGasPressure_, maxGasPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!gasThermalConductivityInitialized_)
                {
                    gasThermalConductivityInitialized_ = tabularizeGasThermalConductivity_();
                    if (gasThermalConductivityInitialized_)
                        return gasThermalConductivity(temperature, pressure);
                }
            }

            printWarning_("gasThermalConductivity", temperature, pressure);
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
        Scalar result = interpolateTP_(liquidThermalConductivity_, temperature, pressure,
                                       pressLiquidIdx_, minLiquidPressure_, maxLiquidPressure_);
        using std::isnan;
        if (isnan(result))
        {
            if constexpr (Multithreading::isSerial())
            {
                // lazy tabularization
                if (!liquidThermalConductivityInitialized_)
                {
                    liquidThermalConductivityInitialized_ = tabularizeLiquidThermalConductivity_();
                    if (liquidThermalConductivityInitialized_)
                        return liquidThermalConductivity(temperature, pressure);
                }
            }

            printWarning_("liquidThermalConductivity", temperature, pressure);
            return RawComponent::liquidThermalConductivity(temperature, pressure);
        }
        return result;
    }


private:
    //! prints a warning if the result is not in range or the table has not been initialized
    static void printWarning_(const std::string& quantity, Scalar arg1, Scalar arg2)
    {
#ifndef NDEBUG
        if (warningPrinted_)
            return;

        if (!initialized_)
            std::cerr << "Warning: tabulated component '" << name()
                      << "' has not been initialized. "
                      << "Call FluidSystem::init() to use the tabulation in order to reduce runtime. \n";
        else
            std::cerr << "Warning: "<<quantity<<"(T="<<arg1<<", p="<<arg2<<") of component '"<<name()
                      << "' is outside tabulation range: ("<< tempMin_<<"<=T<="<<tempMax_<<"), ("
                      << pressMin_<<"<=p<=" <<pressMax_<<"). "
                      << "Forwarded to FluidSystem for direct evaluation of "<<quantity<<". \n";
        warningPrinted_ = true;
#endif
    }

    //! initializes vapor pressure if useVaporPressure = true
    template< bool useVP = useVaporPressure, std::enable_if_t<useVP, int> = 0 >
    static void initVaporPressure_()
    {
        // fill the temperature-pressure arrays
        Dumux::parallelFor(nTemp_, [&](unsigned iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;
            vaporPressure_[iT] = RawComponent::vaporPressure(temperature);
        });
    }

    //! if !useVaporPressure, do nothing here
    template< bool useVP = useVaporPressure, std::enable_if_t<!useVP, int> = 0 >
    static void initVaporPressure_() {}

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
    template<class PropFunc, class MinPFunc, class MaxPFunc>
    static void initTPArray_(PropFunc&& f, MinPFunc&& minP,  MaxPFunc&& maxP, std::vector<typename RawComponent::Scalar>& values)
    {
        Dumux::parallelFor(nTemp_, [&](unsigned iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            Scalar pMax = maxP(iT);
            Scalar pMin = minP(iT);
            for (unsigned iP = 0; iP < nPress_; ++ iP)
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
    template<class RhoFunc, class MinPFunc, class MaxPFunc>
    static void initMinMaxRhoArray_(RhoFunc&& rho,
                                    MinPFunc&& minP,
                                    MaxPFunc&& maxP,
                                    std::vector<typename RawComponent::Scalar>& rhoMin,
                                    std::vector<typename RawComponent::Scalar>& rhoMax)
    {
        Dumux::parallelFor(nTemp_, [&](unsigned iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            rhoMin[iT] = rho(temperature, minP(iT));
            if (iT < nTemp_ - 1)
                rhoMax[iT] = rho(temperature, maxP(iT + 1));
            else
                rhoMax[iT] = rho(temperature, maxP(iT));
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
    static void initPressureArray_(std::vector<typename RawComponent::Scalar>& pressure, PFunc&& p,
                                   const std::vector<typename RawComponent::Scalar>& rhoMin,
                                   const std::vector<typename RawComponent::Scalar>& rhoMax)
    {
        Dumux::parallelFor(nTemp_, [&](unsigned iT)
        {
            Scalar temperature = iT * (tempMax_ - tempMin_)/(nTemp_ - 1) + tempMin_;

            for (unsigned iRho = 0; iRho < nDensity_; ++ iRho)
            {
                Scalar density = Scalar(iRho)/(nDensity_ - 1)
                                 * (rhoMax[iT] - rhoMin[iT])
                                 +  rhoMin[iT];
                pressure[iT + iRho*nTemp_] = p(temperature, density);
            }
        });
    }

    //! returns an interpolated value depending on temperature
    static Scalar interpolateT_(const std::vector<typename RawComponent::Scalar>& values, Scalar T)
    {
        Scalar alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1)
            return std::numeric_limits<Scalar>::quiet_NaN();

        unsigned iT = (unsigned) alphaT;
        alphaT -= iT;

        return values[iT    ]*(1 - alphaT) +
               values[iT + 1]*(    alphaT);
    }

    //! returns an interpolated value depending on temperature and pressure
    template<class GetPIdx, class MinPFunc, class MaxPFunc>
    static Scalar interpolateTP_(const std::vector<typename RawComponent::Scalar>& values, Scalar T, Scalar p,
                                 GetPIdx&& getPIdx, MinPFunc&& minP, MaxPFunc&& maxP)
    {
        Scalar alphaT = tempIdx_(T);
        if (alphaT < 0 || alphaT >= nTemp_ - 1) {
            return std::numeric_limits<Scalar>::quiet_NaN();
        }
        using std::min;
        using std::max;
        unsigned iT = max<int>(0, min<int>(nTemp_ - 2, (int) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = getPIdx(p, iT);
        Scalar alphaP2 = getPIdx(p, iT + 1);

        unsigned iP1 = max<int>(0, min<int>(nPress_ - 2, (int) alphaP1));
        unsigned iP2 = max<int>(0, min<int>(nPress_ - 2, (int) alphaP2));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

#if 0 && !defined NDEBUG
        if(!(0 <= alphaT && alphaT <= 1.0))
            DUNE_THROW(NumericalProblem, "Temperature out of range: "
                       << "T=" << T << " range: [" << tempMin_ << ", " << tempMax_ << "]");
        if(!(0 <= alphaP1 && alphaP1 <= 1.0))
            DUNE_THROW(NumericalProblem, "First pressure out of range: "
                       << "p=" << p << " range: [" << minP(tempIdx_(T)) << ", " << maxP(tempIdx_(T)) << "]");
        if(!(0 <= alphaP2 && alphaP2 <= 1.0))
            DUNE_THROW(NumericalProblem, "Second pressure out of range: "
                       << "p=" << p << " range: [" << minP(tempIdx_(T) + 1) << ", " << maxP(tempIdx_(T) + 1) << "]");
#endif

        return values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
               values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
               values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
               values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }

    //! returns an interpolated value for gas depending on temperature and density
    template<class GetRhoIdx>
    static Scalar interpolateTRho_(const std::vector<typename RawComponent::Scalar>& values, Scalar T, Scalar rho, GetRhoIdx&& rhoIdx)
    {
        using std::min;
        using std::max;
        Scalar alphaT = tempIdx_(T);
        unsigned iT = max<int>(0, min<int>(nTemp_ - 2, (int) alphaT));
        alphaT -= iT;

        Scalar alphaP1 = rhoIdx(rho, iT);
        Scalar alphaP2 = rhoIdx(rho, iT + 1);
        unsigned iP1 = max<int>(0, min<int>(nDensity_ - 2, (int) alphaP1));
        unsigned iP2 = max<int>(0, min<int>(nDensity_ - 2, (int) alphaP2));
        alphaP1 -= iP1;
        alphaP2 -= iP2;

        return values[(iT    ) + (iP1    )*nTemp_]*(1 - alphaT)*(1 - alphaP1) +
               values[(iT    ) + (iP1 + 1)*nTemp_]*(1 - alphaT)*(    alphaP1) +
               values[(iT + 1) + (iP2    )*nTemp_]*(    alphaT)*(1 - alphaP2) +
               values[(iT + 1) + (iP2 + 1)*nTemp_]*(    alphaT)*(    alphaP2);
    }

    //! returns the index of an entry in a temperature field
    static Scalar tempIdx_(Scalar temperature)
    {
        return (nTemp_ - 1)*(temperature - tempMin_)/(tempMax_ - tempMin_);
    }

    //! returns the index of an entry in a pressure field
    static Scalar pressLiquidIdx_(Scalar pressure, unsigned tempIdx)
    {
        Scalar plMin = minLiquidPressure_(tempIdx);
        Scalar plMax = maxLiquidPressure_(tempIdx);

        return (nPress_ - 1)*(pressure - plMin)/(plMax - plMin);
    }

    //! returns the index of an entry in a temperature field
    static Scalar pressGasIdx_(Scalar pressure, unsigned tempIdx)
    {
        Scalar pgMin = minGasPressure_(tempIdx);
        Scalar pgMax = maxGasPressure_(tempIdx);

        return (nPress_ - 1)*(pressure - pgMin)/(pgMax - pgMin);
    }

    //! returns the index of an entry in a density field
    static Scalar densityLiquidIdx_(Scalar density, unsigned tempIdx)
    {
        Scalar densityMin = minLiquidDensity_[tempIdx];
        Scalar densityMax = maxLiquidDensity_[tempIdx];
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    //! returns the index of an entry in a density field
    static Scalar densityGasIdx_(Scalar density, unsigned tempIdx)
    {
        Scalar densityMin = minGasDensity_[tempIdx];
        Scalar densityMax = maxGasDensity_[tempIdx];
        return (nDensity_ - 1) * (density - densityMin)/(densityMax - densityMin);
    }

    //! returns the minimum tabularized liquid pressure at a given temperature index
    static Scalar minLiquidPressure_(int tempIdx)
    {
        using std::max;
        if (!useVaporPressure)
            return pressMin_;
        else
            return max(pressMin_, vaporPressure_[tempIdx] / 1.1);
    }

    //! returns the maximum tabularized liquid pressure at a given temperature index
    static Scalar maxLiquidPressure_(int tempIdx)
    {
        using std::max;
        if (!useVaporPressure)
            return pressMax_;
        else
            return max(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }

    //! returns the minimum tabularized gas pressure at a given temperature index
    static Scalar minGasPressure_(int tempIdx)
    {
        using std::min;
        if (!useVaporPressure)
            return pressMin_;
        else
            return min(pressMin_, vaporPressure_[tempIdx] / 1.1 );
    }

    //! returns the maximum tabularized gas pressure at a given temperature index
    static Scalar maxGasPressure_(int tempIdx)
    {
        using std::min;
        if (!useVaporPressure)
            return pressMax_;
        else
            return min(pressMax_, vaporPressure_[tempIdx] * 1.1);
    }

    template<class RC = RawComponent>
    static bool tabularizeGasEnthalpy_()
    {
        if constexpr (Detail::hasGasEnthalpy<RC>())
        {
            auto gasEnth = [] (auto T, auto p) { return RC::gasEnthalpy(T, p); };
            initTPArray_(gasEnth, minGasPressure_, maxGasPressure_, gasEnthalpy_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeLiquidEnthalpy_()
    {
        if constexpr (Detail::hasLiquidEnthalpy<RC>())
        {
            auto liqEnth = [] (auto T, auto p) { return RC::liquidEnthalpy(T, p); };
            initTPArray_(liqEnth, minLiquidPressure_, maxLiquidPressure_, liquidEnthalpy_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeGasHeatCapacity_()
    {
        if constexpr (Detail::hasGasHeatCapacity<RC>())
        {
            auto gasHC = [] (auto T, auto p) { return RC::gasHeatCapacity(T, p); };
            initTPArray_(gasHC, minGasPressure_, maxGasPressure_, gasHeatCapacity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeLiquidHeatCapacity_()
    {
        if constexpr (Detail::hasLiquidHeatCapacity<RC>())
        {
            auto liqHC = [] (auto T, auto p) { return RC::liquidHeatCapacity(T, p); };
            initTPArray_(liqHC, minLiquidPressure_, maxLiquidPressure_, liquidHeatCapacity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeMinMaxGasDensity_()
    {
        if constexpr (Detail::hasGasDensity<RC>())
        {
            auto gasRho = [] (auto T, auto p) { return RC::gasDensity(T, p); };
            initMinMaxRhoArray_(gasRho, minGasPressure_, maxGasPressure_, minGasDensity_, maxGasDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeMinMaxLiquidDensity_()
    {
        if constexpr (Detail::hasGasEnthalpy<RC>())
        {
            auto liqRho = [] (auto T, auto p) { return RC::liquidDensity(T, p); };
            initMinMaxRhoArray_(liqRho, minLiquidPressure_, maxLiquidPressure_, minLiquidDensity_, maxLiquidDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeGasPressure_()
    {
        // pressure is only defined if the gas is compressible (this is usually the case)
        if constexpr (Detail::hasGasPressure<RC>() && RC::gasIsCompressible())
        {
            auto gasPFunc = [] (auto T, auto rho) { return RC::gasPressure(T, rho); };
            initPressureArray_(gasPressure_, gasPFunc, minGasDensity_, maxGasDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeLiquidPressure_()
    {
        // pressure is only defined if the liquid is compressible (this is often not the case)
        if constexpr (Detail::hasLiquidPressure<RC>() && RC::liquidIsCompressible())
        {
            auto liqPFunc = [] (auto T, auto rho) { return RC::liquidPressure(T, rho); };
            initPressureArray_(liquidPressure_, liqPFunc, minLiquidDensity_, maxLiquidDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeGasDensity_()
    {
        if constexpr (Detail::hasGasDensity<RC>())
        {
            auto gasRho = [] (auto T, auto p) { return RC::gasDensity(T, p); };
            initTPArray_(gasRho, minGasPressure_, maxGasPressure_, gasDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeLiquidDensity_()
    {
        if constexpr (Detail::hasLiquidDensity<RC>())
        {
            // TODO: we could get rid of the lambdas and pass the functor directly. But,
            //       currently Brine is a component (and not a fluid system) expecting a
            //       third argument with a default, which cannot be wrapped in a function pointer.
            //       For this reason we have to wrap this into a lambda here.
            auto liqRho = [] (auto T, auto p) { return RC::liquidDensity(T, p); };
            initTPArray_(liqRho, minLiquidPressure_, maxLiquidPressure_, liquidDensity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeGasViscosity_()
    {
        if constexpr (Detail::hasGasViscosity<RC>())
        {
            auto gasVisc = [] (auto T, auto p) { return RC::gasViscosity(T, p); };
            initTPArray_(gasVisc, minGasPressure_, maxGasPressure_, gasViscosity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeLiquidViscosity_()
    {
        if constexpr (Detail::hasLiquidViscosity<RC>())
        {
            auto liqVisc = [] (auto T, auto p) { return RC::liquidViscosity(T, p); };
            initTPArray_(liqVisc, minLiquidPressure_, maxLiquidPressure_, liquidViscosity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeGasThermalConductivity_()
    {
        if constexpr (Detail::hasGasThermalConductivity<RC>())
        {
            auto gasTC = [] (auto T, auto p) { return RC::gasThermalConductivity(T, p); };
            initTPArray_(gasTC, minGasPressure_, maxGasPressure_, gasThermalConductivity_);
            return true;
        }

        return false;
    }

    template<class RC = RawComponent>
    static bool tabularizeLiquidThermalConductivity_()
    {
        if constexpr (Detail::hasLiquidThermalConductivity<RC>())
        {
            auto liqTC = [] (auto T, auto p) { return RC::liquidThermalConductivity(T, p); };
            initTPArray_(liqTC, minLiquidPressure_, maxLiquidPressure_, liquidThermalConductivity_);
            return true;
        }

        return false;
    }

    // specifies whether the table was initialized
    static bool initialized_;

#ifndef NDEBUG
    // specifies whether some warning was printed
    static bool warningPrinted_;
#endif

    // 1D fields with the temperature as degree of freedom
    static std::vector<typename RawComponent::Scalar> vaporPressure_;

    static std::vector<typename RawComponent::Scalar> minLiquidDensity_;
    static std::vector<typename RawComponent::Scalar> maxLiquidDensity_;
    static bool minMaxLiquidDensityInitialized_;

    static std::vector<typename RawComponent::Scalar> minGasDensity_;
    static std::vector<typename RawComponent::Scalar> maxGasDensity_;
    static bool minMaxGasDensityInitialized_;

    // 2D fields with the temperature and pressure as degrees of freedom
    static std::vector<typename RawComponent::Scalar> gasEnthalpy_;
    static std::vector<typename RawComponent::Scalar> liquidEnthalpy_;
    static bool gasEnthalpyInitialized_;
    static bool liquidEnthalpyInitialized_;

    static std::vector<typename RawComponent::Scalar> gasHeatCapacity_;
    static std::vector<typename RawComponent::Scalar> liquidHeatCapacity_;
    static bool gasHeatCapacityInitialized_;
    static bool liquidHeatCapacityInitialized_;

    static std::vector<typename RawComponent::Scalar> gasDensity_;
    static std::vector<typename RawComponent::Scalar> liquidDensity_;
    static bool gasDensityInitialized_;
    static bool liquidDensityInitialized_;

    static std::vector<typename RawComponent::Scalar> gasViscosity_;
    static std::vector<typename RawComponent::Scalar> liquidViscosity_;
    static bool gasViscosityInitialized_;
    static bool liquidViscosityInitialized_;

    static std::vector<typename RawComponent::Scalar> gasThermalConductivity_;
    static std::vector<typename RawComponent::Scalar> liquidThermalConductivity_;
    static bool gasThermalConductivityInitialized_;
    static bool liquidThermalConductivityInitialized_;

    // 2D fields with the temperature and density as degrees of freedom
    static std::vector<typename RawComponent::Scalar> gasPressure_;
    static std::vector<typename RawComponent::Scalar> liquidPressure_;
    static bool gasPressureInitialized_;
    static bool liquidPressureInitialized_;

    // temperature, pressure and density ranges
    static Scalar tempMin_;
    static Scalar tempMax_;
    static unsigned nTemp_;

    static Scalar pressMin_;
    static Scalar pressMax_;
    static unsigned nPress_;

    static Scalar densityMin_;
    static Scalar densityMax_;
    static unsigned nDensity_;
};


template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::initialized_ = false;

#ifndef NDEBUG
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::warningPrinted_ = false;
#endif

template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::minMaxLiquidDensityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::minMaxGasDensityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::gasEnthalpyInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::liquidEnthalpyInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::gasHeatCapacityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::liquidHeatCapacityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::gasDensityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::liquidDensityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::gasViscosityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::liquidViscosityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::gasThermalConductivityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::liquidThermalConductivityInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::gasPressureInitialized_ = false;
template <class RawComponent, bool useVaporPressure>
bool TabulatedComponent<RawComponent, useVaporPressure>::liquidPressureInitialized_ = false;

template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::vaporPressure_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::minLiquidDensity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::maxLiquidDensity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::minGasDensity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::maxGasDensity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::gasEnthalpy_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::liquidEnthalpy_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::gasHeatCapacity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::liquidHeatCapacity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::gasDensity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::liquidDensity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::gasViscosity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::liquidViscosity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::gasThermalConductivity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::liquidThermalConductivity_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::gasPressure_;
template <class RawComponent, bool useVaporPressure>
std::vector<typename RawComponent::Scalar> TabulatedComponent<RawComponent, useVaporPressure>::liquidPressure_;
template <class RawComponent, bool useVaporPressure>
typename RawComponent::Scalar TabulatedComponent<RawComponent, useVaporPressure>::tempMin_;
template <class RawComponent, bool useVaporPressure>
typename RawComponent::Scalar TabulatedComponent<RawComponent, useVaporPressure>::tempMax_;
template <class RawComponent, bool useVaporPressure>
unsigned TabulatedComponent<RawComponent, useVaporPressure>::nTemp_;
template <class RawComponent, bool useVaporPressure>
typename RawComponent::Scalar TabulatedComponent<RawComponent, useVaporPressure>::pressMin_;
template <class RawComponent, bool useVaporPressure>
typename RawComponent::Scalar TabulatedComponent<RawComponent, useVaporPressure>::pressMax_;
template <class RawComponent, bool useVaporPressure>
unsigned TabulatedComponent<RawComponent, useVaporPressure>::nPress_;
template <class RawComponent, bool useVaporPressure>
typename RawComponent::Scalar TabulatedComponent<RawComponent, useVaporPressure>::densityMin_;
template <class RawComponent, bool useVaporPressure>
typename RawComponent::Scalar TabulatedComponent<RawComponent, useVaporPressure>::densityMax_;
template <class RawComponent, bool useVaporPressure>
unsigned TabulatedComponent<RawComponent, useVaporPressure>::nDensity_;

// forward declaration
template <class Component>
struct IsAqueous;

// we are aqueous if the raw compont is so
template <class RawComponent, bool useVaporPressure>
struct IsAqueous<TabulatedComponent<RawComponent, useVaporPressure>> : public IsAqueous<RawComponent> {};

} // end namespace Dumux::Components

#endif
