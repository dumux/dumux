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
 * \ingroup MaterialTests
 * \brief Plot properties of components and fluids.
 */

#include "config.h"
#include <array>
#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>
#include <dumux/common/typetraits/isvalid.hh>
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
struct checkLiqDen { template<class C> auto operator()(C&& c) -> decltype(C::template liquidDensity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkLiqEnth { template<class C> auto operator()(C&& c) -> decltype(C::template liquidEnthalpy<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkLiqHeatCap { template<class C> auto operator()(C&& c) -> decltype(C::template liquidHeatCapacity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkLiqVisc { template<class C> auto operator()(C&& c) -> decltype(C::template liquidViscosity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkLiqThermCond { template<class C> auto operator()(C&& c) -> decltype(C::template liquidThermalConductivity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkGasDen { template<class C> auto operator()(C&& c) -> decltype(C::template gasDensity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkGasEnth { template<class C> auto operator()(C&& c) -> decltype(C::template gasEnthalpy<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkGasHeatCap { template<class C> auto operator()(C&& c) -> decltype(C::template gasHeatCapacity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkGasVisc { template<class C> auto operator()(C&& c) -> decltype(C::template gasViscosity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkGasThermCond { template<class C> auto operator()(C&& c) -> decltype(C::template gasThermalConductivity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkSolDen { template<class C> auto operator()(C&& c) -> decltype(C::template solidDensity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkSolHeatCap { template<class C> auto operator()(C&& c) -> decltype(C::template solidHeatCapacity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkSolThermCond { template<class C> auto operator()(C&& c) -> decltype(C::template solidThermalConductivity<DisableStaticAssert>(0.0, 0.0)) {} };
struct checkIonCharge { template<class C> auto operator()(C&& c) -> decltype(C::template charge<DisableStaticAssert>(0.0, 0.0)) {} };

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

//! Plot properties if overloads compile
template<class C, class hasNoDensityOverload = checkLiqDen>
auto plotLiquidDensity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoDensityOverload{})(declval<C>()))::value && ComponentTraits<C>::hasLiquidState, void>
{
    auto f = [] (auto T, auto p) { return C::liquidDensity(T, p); };
    plot(f, T, p, C::name(), "liquid", "density", "[kg/^3]", openPlot);
}

template<class C, class hasNoEnthalpyOverload = checkLiqEnth>
auto plotLiquidEnthalpy(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoEnthalpyOverload{})(declval<C>()))::value && ComponentTraits<C>::hasLiquidState, void>
{
    auto f = [] (auto T, auto p) { return C::liquidEnthalpy(T, p); };
    plot(f, T, p, C::name(), "liquid", "enthalpy", "[J/(kg)]", openPlot);
}

template<class C, class hasNoHeatCapOverload = checkLiqHeatCap>
auto plotLiquidHeatCapacity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoHeatCapOverload{})(declval<C>()))::value && ComponentTraits<C>::hasLiquidState, void>
{
    auto f = [] (auto T, auto p) { return C::liquidHeatCapacity(T, p); };
    plot(f, T, p, C::name(), "liquid", "heat capacity", "[J/(kg*K)]", openPlot);
}

template<class C, class hasNoViscOverload = checkLiqVisc>
auto plotLiquidViscosity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoViscOverload{})(declval<C>()))::value && ComponentTraits<C>::hasLiquidState, void>
{
    auto f = [] (auto T, auto p) { return C::liquidViscosity(T, p); };
    plot(f, T, p, C::name(), "liquid", "viscosity", "[Pa*s]", openPlot);
}

template<class C, class hasNoThermCondOverload = checkLiqThermCond>
auto plotLiquidThermalConductivity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoThermCondOverload{})(declval<C>()))::value && ComponentTraits<C>::hasLiquidState, void>
{
    auto f = [] (auto T, auto p) { return C::liquidThermalConductivity(T, p); };
    plot(f, T, p, C::name(), "liquid", "thermal conductivity", "[J/(kg*K)]", openPlot);
}

template<class C, class hasNoDensityOverload = checkGasDen>
auto plotGasDensity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoDensityOverload{})(declval<C>()))::value && ComponentTraits<C>::hasGasState, void>
{
    auto f = [] (auto T, auto p) { return C::gasDensity(T, p); };
    plot(f, T, p, C::name(), "gas", "density", "[kg/^3]", openPlot);
}

template<class C, class hasNoEnthalpyOverload = checkGasEnth>
auto plotGasEnthalpy(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoEnthalpyOverload{})(declval<C>()))::value && ComponentTraits<C>::hasGasState, void>
{
    auto f = [] (auto T, auto p) { return C::gasEnthalpy(T, p); };
    plot(f, T, p, C::name(), "gas", "enthalpy", "[J/(kg)]", openPlot);
}

template<class C, class hasNoHeatCapOverload = checkGasHeatCap>
auto plotGasHeatCapacity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoHeatCapOverload{})(declval<C>()))::value && ComponentTraits<C>::hasGasState, void>
{
    auto f = [] (auto T, auto p) { return C::gasHeatCapacity(T, p); };
    plot(f, T, p, C::name(), "gas", "heat capacity", "[J/(kg*K)]", openPlot);
}

template<class C, class hasNoViscOverload = checkGasVisc>
auto plotGasViscosity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoViscOverload{})(declval<C>()))::value && ComponentTraits<C>::hasGasState, void>
{
    auto f = [] (auto T, auto p) { return C::gasViscosity(T, p); };
    plot(f, T, p, C::name(), "gas", "viscosity", "[Pa*s]", openPlot);
}

template<class C, class hasNoThermCondOverload = checkGasThermCond>
auto plotGasThermalConductivity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoThermCondOverload{})(declval<C>()))::value && ComponentTraits<C>::hasGasState, void>
{
    auto f = [] (auto T, auto p) { return C::gasThermalConductivity(T, p); };
    plot(f, T, p, C::name(), "gas", "thermal conductivity", "[J/(kg*K)]", openPlot);
}

template<class C, class hasNoDensityOverload = checkSolDen>
auto plotSolidDensity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoDensityOverload{})(declval<C>()))::value && ComponentTraits<C>::hasSolidState, void>
{
    auto f = [] (auto T, auto p) { return C::solidDensity(T); };
    plot(f, T, p, C::name(), "solid", "density", "[kg/^3]", openPlot);
}

template<class C, class hasNoHeatCapOverload = checkSolHeatCap>
auto plotSolidHeatCapacity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoHeatCapOverload{})(declval<C>()))::value && ComponentTraits<C>::hasSolidState, void>
{
    auto f = [] (auto T, auto p) { return C::solidHeatCapacity(T); };
    plot(f, T, p, C::name(), "solid", "heat capacity", "[J/(kg*K)]", openPlot);
}

template<class C, class hasNoThermCondOverload = checkSolThermCond>
auto plotSolidThermalConductivity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoThermCondOverload{})(declval<C>()))::value && ComponentTraits<C>::hasSolidState, void>
{
    auto f = [] (auto T, auto p) { return C::solidThermalConductivity(T); };
    plot(f, T, p, C::name(), "solid", "thermal conductivity", "[J/(kg*K)]", openPlot);
}

template<class C, class hasNoChargeOverload = checkIonCharge>
auto plotIonCharge(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<!decltype(isValid(hasNoChargeOverload{})(declval<C>()))::value && ComponentTraits<C>::isIon, void>
{
    auto f = [] (auto T, auto p) { return C::charge(); };
    plot(f, T, p, C::name(), "ion", "charge", "[e]", openPlot);
}

//! Do not plot properties if overloads don't compile
template<class C, class hasNoDensityOverload = checkLiqDen>
auto plotLiquidDensity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoDensityOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasLiquidState, void> {}

template<class C, class hasNoEnthalpyOverload = checkLiqEnth>
auto plotLiquidEnthalpy(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoEnthalpyOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasLiquidState, void> {}

template<class C, class hasNoHeatCapOverload = checkLiqHeatCap>
auto plotLiquidHeatCapacity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoHeatCapOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasLiquidState, void> {}

template<class C, class hasNoViscOverload = checkLiqVisc>
auto plotLiquidViscosity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoViscOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasLiquidState, void> {}

template<class C, class hasNoThermCondOverload = checkLiqThermCond>
auto plotLiquidThermalConductivity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoThermCondOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasLiquidState, void> {}

template<class C, class hasNoDensityOverload = checkGasDen>
auto plotGasDensity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoDensityOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasGasState, void> {}

template<class C, class hasNoEnthalpyOverload = checkGasEnth>
auto plotGasEnthalpy(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoEnthalpyOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasGasState, void> {}

template<class C, class hasNoHeatCapOverload = checkGasHeatCap>
auto plotGasHeatCapacity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoHeatCapOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasGasState, void> {}

template<class C, class hasNoViscOverload = checkGasVisc>
auto plotGasViscosity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoViscOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasGasState, void> {}

template<class C, class hasNoThermCondOverload = checkGasThermCond>
auto plotGasThermalConductivity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoThermCondOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasGasState, void> {}

template<class C, class hasNoDensityOverload = checkSolDen>
auto plotSolidDensity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoDensityOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasSolidState, void> {}

template<class C, class hasNoHeatCapOverload = checkSolHeatCap>
auto plotSolidHeatCapacity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoHeatCapOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasSolidState, void> {}

template<class C, class hasNoThermCondOverload = checkSolThermCond>
auto plotSolidThermalConductivity(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoThermCondOverload{})(declval<C>()))::value || !ComponentTraits<C>::hasSolidState, void> {}

template<class C, class hasNoChargeOverload = checkIonCharge>
auto plotIonCharge(const vector<double>& T, double p, bool openPlot)
-> typename std::enable_if_t<decltype(isValid(hasNoChargeOverload{})(declval<C>()))::value || !ComponentTraits<C>::isIon, void> {}

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
