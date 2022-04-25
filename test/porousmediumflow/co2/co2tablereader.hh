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
 * \brief A generic template for tabulated material laws that depend
 *        on two parameters.
 */
#ifndef DUMUX_TABULATED_CO2_HH
#define DUMUX_TABULATED_CO2_HH

#warning "This header is deprecated and will be removed after release 3.5. Copy co2tables.hh and maybe adjust include 'inc' include to your needs."
#include <dune/common/float_cmp.hh>
#include <dumux/common/exceptions.hh>

namespace Dumux {
/*!
 * \ingroup Components
 * \brief A generic template for tabulated material laws that depend
 *        on two parameters.
 */
template <class Traits>
class TabulatedCO2Properties
{
    using Scalar = typename Traits::Scalar;
    enum { numTempSteps = Traits::numTempSteps, numPressSteps = Traits::numPressSteps };

public:
    TabulatedCO2Properties() = default;

    constexpr Scalar minTemp() const
    { return Traits::minTemp; }

    constexpr Scalar maxTemp() const
    { return Traits::maxTemp; }

    constexpr Scalar minPress() const
    { return Traits::minPress; }

    constexpr Scalar maxPress() const
    { return Traits::maxPress; }

    constexpr bool applies(Scalar temperature, Scalar pressure) const
    {
        return minTemp() <= temperature && temperature <= maxTemp() &&
            minPress() <= pressure && pressure <= maxPress();
    }

    constexpr Scalar at(Scalar temperature, Scalar pressure) const
    {
        if (!applies(temperature,pressure))
        {
            if (temperature<minTemp())
                temperature=minTemp();
            if(temperature>maxTemp())
                temperature=maxTemp();
            if(pressure<minPress())
                pressure=minPress();
            if(pressure>maxPress())
                pressure=maxPress();
        }

        int i = findTempIdx_(temperature);
        int j = findPressIdx_(pressure);

        Scalar tempAtI = temperatureAt_(i);
        Scalar tempAtI1 = temperatureAt_(i + 1);
        Scalar pressAtI = pressureAt_(j);
        Scalar pressAtI1 = pressureAt_(j + 1);

        Scalar alpha = (temperature - tempAtI)/(tempAtI1 - tempAtI);
        Scalar beta = (pressure - pressAtI)/(pressAtI1 - pressAtI);

        // bi-linear interpolation
        Scalar lowresValue =
            (1-alpha)*(1-beta)*val(i, j) +
            (1-alpha)*(  beta)*val(i, j + 1) +
            (  alpha)*(1-beta)*val(i + 1, j) +
            (  alpha)*(  beta)*val(i + 1, j + 1);

        // return the weighted sum of the low- and high-resolution
        // values
        return lowresValue;
    }

    constexpr Scalar val(int i, int j) const
    {
        return Traits::vals[i][j];
    }

protected:
    constexpr int findTempIdx_(Scalar temperature) const
    {
        if (Dune::FloatCmp::eq<Scalar>(temperature, maxTemp()))
            return numTempSteps - 2;
        const int result = static_cast<int>((temperature - minTemp())/(maxTemp() - minTemp())*(numTempSteps - 1));

        using std::min;
        using std::max;
        return max(0, min(result, numTempSteps - 2));
    }

    constexpr int findPressIdx_(Scalar pressure) const
    {
        if (Dune::FloatCmp::eq<Scalar>(pressure, maxPress()))
            return numPressSteps - 2;
        const int result = static_cast<int>((pressure - minPress())/(maxPress() - minPress())*(numPressSteps - 1));

        using std::min;
        using std::max;
        return max(0, min(result, numPressSteps - 2));
    }

    constexpr Scalar temperatureAt_(int i) const
    { return i*(maxTemp() - minTemp())/(numTempSteps - 1) + minTemp(); }
    constexpr Scalar pressureAt_(int j) const
    { return j*(maxPress() - minPress())/(numPressSteps - 1) + minPress(); }
};

} // end namespace Dumux

#endif
