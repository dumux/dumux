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
 * \brief A reader and tables for CO$_2$ tabulated material laws that depend
 *        on pressure and temperature.
 */
#ifndef DUMUX_EXAMPLE_BIOMINERALIZATION_CO2TABLES_HH
#define DUMUX_EXAMPLE_BIOMINERALIZATION_CO2TABLES_HH

#include <dune/common/float_cmp.hh>

// ## The CO2 tables (`co2tables.hh`)
//
// This file contains the __co2table class__ which forwards to tabulated properties of CO2.
// The tables are generated using the NIST (National Institute of Standards
// and Technology) Standard Reference Database Number 69
// (https://doi.org/10.18434/T4D303).
//
// Copyright for NIST Standard Reference Data is governed by the Standard
// Reference Data Act (https://www.nist.gov/srd/public-law).
//
// ######################################################################
// In case you are using this the data generated with this script
// please cite the following publications:
//
// P.J. Linstrom and W.G. Mallard, Eds.,
// NIST Chemistry WebBook, NIST Standard Reference Database Number 69,
// National Institute of Standards and Technology, Gaithersburg MD, 20899,
// https://doi.org/10.18434/T4D303, (retrieved [insert date]).
//
// Span, Roland, and Wolfgang Wagner.
// "A new equation of state for carbon dioxide covering
// the fluid region from the triple‐point temperature
// to 1100 K at pressures up to 800 MPa."
// Journal of physical and chemical reference data 25.6 (1996): 1509-1596.
// https://doi.org/10.1063/1.555991
//
// ######################################################################
//
// The density and the enthalpy are calculated using the equation of Span and
// Wagner (2009 "A New Equation of State for Carbon Dioxide Covering the Fluid
// Region from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa").
// Therefore, the maximum pressure limit is the lowest of the following values:
// * 800.0000 MPa
// * The pressure at which a density of 1178.5 kg/m3 is reached.
//
// [[content]]
//
// [[codeblock]]

namespace Dumux::BiomineralizationCO2Tables {

/*!
 * \ingroup Components
 * \brief A generic template for tabulated material laws that depend
 *        on two parameters.
 */
template <class Traits>
class TabulatedProperties
{
    using Scalar = typename Traits::Scalar;

    static constexpr auto numTempSteps = Traits::numTempSteps;
    static constexpr auto numPressSteps = Traits::numPressSteps;

public:
    TabulatedProperties() = default;

    constexpr Scalar minTemp() const { return Traits::minTemp; }
    constexpr Scalar maxTemp() const { return Traits::maxTemp; }
    constexpr Scalar minPress() const { return Traits::minPress; }
    constexpr Scalar maxPress() const { return Traits::maxPress; }

    constexpr bool applies(Scalar temperature, Scalar pressure) const
    {
        return minTemp() <= temperature && temperature <= maxTemp() &&
               minPress() <= pressure && pressure <= maxPress();
    }

    constexpr Scalar at(Scalar temperature, Scalar pressure) const
    {
        if (!applies(temperature, pressure))
        {
            if (temperature<minTemp()) temperature = minTemp();
            else if (temperature>maxTemp()) temperature = maxTemp();

            if (pressure<minPress()) pressure = minPress();
            else if (pressure>maxPress()) pressure = maxPress();
        }

        const int i = findTempIdx_(temperature);
        const int j = findPressIdx_(pressure);

        const Scalar tempAtI = temperatureAt_(i);
        const Scalar tempAtI1 = temperatureAt_(i + 1);
        const Scalar pressAtI = pressureAt_(j);
        const Scalar pressAtI1 = pressureAt_(j + 1);

        const Scalar alpha = (temperature - tempAtI)/(tempAtI1 - tempAtI);
        const Scalar beta = (pressure - pressAtI)/(pressAtI1 - pressAtI);

        // bi-linear interpolation
        const Scalar lowresValue =
            (1-alpha)*(1-beta)*val(i, j) +
            (1-alpha)*(  beta)*val(i, j + 1) +
            (  alpha)*(1-beta)*val(i + 1, j) +
            (  alpha)*(  beta)*val(i + 1, j + 1);

        // return the weighted sum of the low- and high-resolution values
        return lowresValue;
    }

    constexpr Scalar val(int i, int j) const
    { return Traits::vals[i][j]; }

private:
    constexpr int findTempIdx_(Scalar temperature) const
    {
        if (Dune::FloatCmp::eq<Scalar>(temperature, maxTemp()))
            return numTempSteps - 2;

        const int result = static_cast<int>((temperature - minTemp())/(maxTemp() - minTemp())*(numTempSteps - 1));

        using std::clamp;
        return clamp(result, 0, numTempSteps - 2);
    }

    constexpr int findPressIdx_(Scalar pressure) const
    {
        if (Dune::FloatCmp::eq<Scalar>(pressure, maxPress()))
            return numPressSteps - 2;

        const int result = static_cast<int>((pressure - minPress())/(maxPress() - minPress())*(numPressSteps - 1));

        using std::clamp;
        return clamp(result, 0, numPressSteps - 2);
    }

    constexpr Scalar temperatureAt_(int i) const
    { return i*(maxTemp() - minTemp())/(numTempSteps - 1) + minTemp(); }

    constexpr Scalar pressureAt_(int j) const
    { return j*(maxPress() - minPress())/(numPressSteps - 1) + minPress(); }
};

#ifndef DOXYGEN // hide from doxygen
// the real work is done by some external program which provides
// ready-to-use tables.
#include "co2values.inc"
#endif

using TabulatedDensity = TabulatedProperties<TabulatedDensityTraits>;
using TabulatedEnthalpy = TabulatedProperties<TabulatedEnthalpyTraits>;

// this class collects all the tabulated quantities in one convenient place
struct CO2Tables
{
   static constexpr inline TabulatedEnthalpy tabulatedEnthalpy = {};
   static constexpr inline TabulatedDensity tabulatedDensity = {};
};

} // end namespace Dumux::GeneratedCO2Tables
// [[/codeblock]]
// [[/content]]

#endif
