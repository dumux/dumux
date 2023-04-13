// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief Plot properties of components and fluids.
 */

#include "config.h"
#include <array>
#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

#include <dune/common/std/type_traits.hh>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/ammonia.hh>
#include <dumux/material/components/benzene.hh>
#include <dumux/material/components/brine.hh>
#include <dumux/material/components/calcite.hh>
#include <dumux/material/components/calciumion.hh>
#include <dumux/material/components/cao.hh>
#include <dumux/material/components/cao2h2.hh>
#include <dumux/material/components/carbonateion.hh>
#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/chlorideion.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/glucose.hh>
#include <dumux/material/components/granite.hh>
#include <dumux/material/components/h2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/heavyoil.hh>
#include <dumux/material/components/mesitylene.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/simpleco2.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/sodiumion.hh>
#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/urea.hh>
#include <dumux/material/components/xylene.hh>
#include <dumux/material/components/componenttraits.hh>

using namespace std;
using namespace Dumux;

namespace Dumux {

//! Helper struct to deactivate static assertions in component's base classes.
struct DisableStaticAssert {};

/*!
 * \brief Specialization of Dumux::AlwaysFalse for the struct defined
 *        above. This is done in order to deactivate the static_assert in
 *        the base classes of components. If the base class function is compiled
 *        we do not call it (see below).
 */
template<>
struct AlwaysFalse<DisableStaticAssert> : public std::true_type {};
} // end namespace Dumux

//! Helper structs for detecting if a component has certain functions overloaded
template<class C> using DetectLiqDen = decltype(C::template liquidDensity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectLiqEnth = decltype(C::template liquidEnthalpy<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectHeatCap = decltype(C::template liquidHeatCapacity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectLiqVisc = decltype(C::template liquidViscosity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectLiqThermCond = decltype(C::template liquidThermalConductivity<DisableStaticAssert>(0.0, 0.0));

template<class C> using DetectGasDen = decltype(C::template gasDensity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectGasEnth = decltype(C::template gasEnthalpy<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectGasHeatCap = decltype(C::template gasHeatCapacity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectGasVisc = decltype(C::template gasViscosity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectGasThermCond = decltype(C::template gasThermalConductivity<DisableStaticAssert>(0.0, 0.0));

template<class C> using DetectSolDen = decltype(C::template solidDensity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectSolHeatCap = decltype(C::template solidHeatCapacity<DisableStaticAssert>(0.0, 0.0));
template<class C> using DetectSolThermCond = decltype(C::template solidThermalConductivity<DisableStaticAssert>(0.0, 0.0));

template<class C> using DetectIonCharge = decltype(C::template charge<DisableStaticAssert>(0.0, 0.0));

//! Plot given values
template<class Functor>
void plot(Functor&& f,
          const vector<double>& T,
          const double pressure,
          const std::string& compName,
          const std::string& phaseName,
          const std::string& propName,
          const std::string& unit,
          bool openPlot)
{
    vector<double> values(T.size());
    for (int i = 0; i < T.size(); ++i)
        values[i] = f(T[i], pressure);

    const auto minMax = minmax_element(values.begin(), values.end());
    Dumux::GnuplotInterface<double> gnuplot(true);
    gnuplot.setOpenPlotWindow(openPlot);
    gnuplot.setCreateImage(true);
    gnuplot.setXRange(T[0], T[T.size()-1]);
    gnuplot.setYRange(*(minMax.first)*0.999, *(minMax.second)*1.001);
    gnuplot.setXlabel("temperature [K]");
    gnuplot.setYlabel(phaseName + " " + propName + " " + unit);
    gnuplot.setDatafileSeparator(',');
    gnuplot.addDataSetToPlot(T, values, compName + "_" + phaseName + "_" + propName + ".csv");
    gnuplot.plot(compName + "_" + phaseName + "_" + propName);
}

template<class C>
void plotLiquidDensity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasLiquidState && !Dune::Std::is_detected<DetectLiqDen, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::liquidDensity(T, p); };
        plot(f, T, p, C::name(), "liquid", "density", "[kg/^3]", openPlot);
    }
}

template<class C>
void plotLiquidEnthalpy(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasLiquidState && !Dune::Std::is_detected<DetectLiqEnth, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::liquidEnthalpy(T, p); };
        plot(f, T, p, C::name(), "liquid", "enthalpy", "[J/(kg)]", openPlot);
    }
}

template<class C>
void plotLiquidHeatCapacity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasLiquidState && !Dune::Std::is_detected<DetectHeatCap, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::liquidHeatCapacity(T, p); };
        plot(f, T, p, C::name(), "liquid", "heat capacity", "[J/(kg*K)]", openPlot);
    }
}

template<class C>
void plotLiquidViscosity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasLiquidState && !Dune::Std::is_detected<DetectLiqVisc, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::liquidViscosity(T, p); };
        plot(f, T, p, C::name(), "liquid", "viscosity", "[Pa*s]", openPlot);
    }
}

template<class C>
void plotLiquidThermalConductivity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasLiquidState && !Dune::Std::is_detected<DetectLiqThermCond, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::liquidThermalConductivity(T, p); };
        plot(f, T, p, C::name(), "liquid", "thermal conductivity", "[J/(kg*K)]", openPlot);
    }
}

template<class C>
void plotGasDensity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasGasState && !Dune::Std::is_detected<DetectGasDen, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::gasDensity(T, p); };
        plot(f, T, p, C::name(), "gas", "density", "[kg/^3]", openPlot);
    }
}

template<class C>
void plotGasEnthalpy(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasGasState && !Dune::Std::is_detected<DetectGasEnth, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::gasEnthalpy(T, p); };
        plot(f, T, p, C::name(), "gas", "enthalpy", "[J/(kg)]", openPlot);
    }
}

template<class C>
void plotGasHeatCapacity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasGasState && !Dune::Std::is_detected<DetectGasHeatCap, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::gasHeatCapacity(T, p); };
        plot(f, T, p, C::name(), "gas", "heat capacity", "[J/(kg*K)]", openPlot);
    }
}

template<class C>
void plotGasViscosity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasGasState && !Dune::Std::is_detected<DetectGasVisc, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::gasViscosity(T, p); };
        plot(f, T, p, C::name(), "gas", "viscosity", "[Pa*s]", openPlot);
    }
}

template<class C>
void plotGasThermalConductivity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasGasState && !Dune::Std::is_detected<DetectGasThermCond, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::gasThermalConductivity(T, p); };
        plot(f, T, p, C::name(), "gas", "thermal conductivity", "[J/(kg*K)]", openPlot);
    }
}

template<class C>
void plotSolidDensity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasSolidState && !Dune::Std::is_detected<DetectSolDen, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::solidDensity(T); };
        plot(f, T, p, C::name(), "solid", "density", "[kg/^3]", openPlot);
    }
}

template<class C>
void plotSolidHeatCapacity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasSolidState && !Dune::Std::is_detected<DetectSolHeatCap, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::solidHeatCapacity(T); };
        plot(f, T, p, C::name(), "solid", "heat capacity", "[J/(kg*K)]", openPlot);
    }
}

template<class C>
void plotSolidThermalConductivity(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::hasSolidState && !Dune::Std::is_detected<DetectSolThermCond, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::solidThermalConductivity(T); };
        plot(f, T, p, C::name(), "solid", "thermal conductivity", "[J/(kg*K)]", openPlot);
    }
}

template<class C>
void plotIonCharge(const vector<double>& T, double p, bool openPlot)
{
    if constexpr (ComponentTraits<C>::isIon && !Dune::Std::is_detected<DetectIonCharge, C>::value)
    {
        auto f = [] (auto T, auto p) { return C::charge(); };
        plot(f, T, p, C::name(), "ion", "charge", "[e]", openPlot);
    }
}

//! A number of properties of a component
template<class Component>
void plotStuff(bool openPlotWindow)
{
    double pressure = 1e5;
    double TMin = 273.15;
    double TMax = 323.15;
    double TRange = TMax - TMin;
    const unsigned int numIntervals = 100;
    vector<double> T(numIntervals + 1);
    for (int i = 0; i <= numIntervals; i++)
        T[i] = TMin + TRange * double(i) /double(numIntervals);

    plotLiquidDensity<Component>(T, pressure, openPlotWindow);
    plotLiquidEnthalpy<Component>(T, pressure, openPlotWindow);
    plotLiquidHeatCapacity<Component>(T, pressure, openPlotWindow);
    plotLiquidViscosity<Component>(T, pressure, openPlotWindow);
    plotLiquidThermalConductivity<Component>(T, pressure, openPlotWindow);

    plotGasDensity<Component>(T, pressure, openPlotWindow);
    plotGasEnthalpy<Component>(T, pressure, openPlotWindow);
    plotGasHeatCapacity<Component>(T, pressure, openPlotWindow);
    plotGasViscosity<Component>(T, pressure, openPlotWindow);
    plotGasThermalConductivity<Component>(T, pressure, openPlotWindow);

    plotSolidDensity<Component>(T, pressure, openPlotWindow);
    plotSolidThermalConductivity<Component>(T, pressure, openPlotWindow);
    plotSolidHeatCapacity<Component>(T, pressure, openPlotWindow);

    plotIonCharge<Component>(T, pressure, openPlotWindow);
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char *argv[])
{
    using namespace Dumux;

    bool openPlotWindow = false;
    if (argc == 3 && (strcmp(argv[2], "1") || strcmp(argv[2], "true") || strcmp(argv[2], "True")))
        openPlotWindow = true;

    if (argc < 2)
        DUNE_THROW(Dune::InvalidStateException, "At least one argument (the component name) is required!");

    const std::string compName = argv[1];

    if (compName == "Air")
        plotStuff< Components::Air<double> >(openPlotWindow);
    else if (compName == "Ammonia")
        plotStuff< Components::Ammonia<double> >(openPlotWindow);
    else if (compName == "Benzene")
        plotStuff< Components::Benzene<double> >(openPlotWindow);
    else if (compName == "Brine")
    {
        Parameters::init([](auto& params){ params["Brine.Salinity"] = "0.1"; });
        plotStuff< Components::Brine<double> >(openPlotWindow);
    }
    else if (compName == "Calcite")
        plotStuff< Components::Calcite<double> >(openPlotWindow);
    else if (compName == "CalciumIon")
        plotStuff< Components::CalciumIon<double> >(openPlotWindow);
    else if (compName == "CaO")
        plotStuff< Components::CaO<double> >(openPlotWindow);
    else if (compName == "CaO2H2")
        plotStuff< Components::CaO2H2<double> >(openPlotWindow);
    else if (compName == "CarbonateIon")
        plotStuff< Components::CarbonateIon<double> >(openPlotWindow);
    else if (compName == "CH4")
        plotStuff< Components::CH4<double> >(openPlotWindow);
    else if (compName == "ChlorideIon")
        plotStuff< Components::ChlorideIon<double> >(openPlotWindow);
    else if (compName == "Constant")
    {
        Parameters::init([](auto& params){
            params["Component.LiquidDensity"] = "1e3";
            params["Component.LiquidKinematicViscosity"] = "1e-3";
            params["Component.LiquidThermalConductivity"] = "0.679";
            params["Component.LiquidHeatCapacity"] = "4.2e3";
            params["Component.GasDensity"] = "1";
            params["Component.GasKinematicViscosity"] = "1";
            params["Component.SolidDensity"] = "1e3";
            params["Component.SolidThermalConductivity"] = "0.679";
            params["Component.SolidHeatCapacity"] = "4.2e3";
            params["Component.EnthalpyOfVaporization"] ="2453e3";
            params["Component.GasThermalConductivity"] ="1";
            params["Component.GasHeatCapacity"] ="1";
        });
        plotStuff< Components::Constant<1, double> >(openPlotWindow);
    }
    else if (compName == "Glucose")
        plotStuff< Components::Glucose<double> >(openPlotWindow);
    else if (compName == "Granite")
        plotStuff< Components::Granite<double> >(openPlotWindow);
    else if (compName == "H2")
        plotStuff< Components::H2<double> >(openPlotWindow);
    else if (compName == "H2O")
        plotStuff< Components::H2O<double> >(openPlotWindow);
    else if (compName == "HeavyOil")
        plotStuff< Components::HeavyOil<double> >(openPlotWindow);
    else if (compName == "Mesitylene")
        plotStuff< Components::Mesitylene<double> >(openPlotWindow);
    else if (compName == "N2")
        plotStuff< Components::N2<double> >(openPlotWindow);
    else if (compName == "NaCl")
        plotStuff< Components::NaCl<double> >(openPlotWindow);
    else if (compName == "O2")
        plotStuff< Components::O2<double> >(openPlotWindow);
    else if (compName == "SimpleCO2")
        plotStuff< Components::SimpleCO2<double>  >(openPlotWindow);
    else if (compName == "SimpleH2O")
        plotStuff< Components::SimpleH2O<double>  >(openPlotWindow);
    else if (compName == "SodiumIon")
        plotStuff< Components::SodiumIon<double>  >(openPlotWindow);
    else if (compName == "Trichloroethene")
        plotStuff< Components::Trichloroethene<double> >(openPlotWindow);
    else if (compName == "Urea")
        plotStuff< Components::Urea<double>  >(openPlotWindow);
    else if (compName == "Xylene")
        plotStuff< Components::Xylene<double> >(openPlotWindow);
    else
        DUNE_THROW(Dune::NotImplemented, "Test for component " << compName);
}
