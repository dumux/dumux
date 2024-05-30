// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Shomate equations for enthalpy and heat capacity.
 *
 * The Shomate equations were inroduced in \cite Shomate1954.
 * They offer analytical equations for computing the heat capacity as well as enthalpy.
 * Both equations use a set of component-specific parameters \f$A,B,C,D,E,F,G,H\f$.
 * For the heat capacity \f$C_p^0\f$, one obtains
 \f[
 C_p^0 = A + Bt + Ct^2 + Dt^3 + E/t^2,
 \f]
 * while for the enthalpy with a reference state at \f$T=298.15\text{K}\f$, one uses
 \f[
 H^0-H_{298.15}^0 = At + \frac{Bt^2}{2} + \frac{Ct^3}{3} + \frac{Dt^4}{4} -\frac{E}{t} + F - H
 \f]
 *where:
 * * \f$ C_p \f$ is the heat capacity in J/(mol*K),
 * * \f$ H^0 \f$ represents the standard enthalpy in kJ/mol,
 * * \f$ t \f$ is the temperature in K divided by 1000.
 *
 */
#ifndef DUMUX_MATERIAL_COMPONENTS_SHOMATE_HH
#define DUMUX_MATERIAL_COMPONENTS_SHOMATE_HH

#include <algorithm>
#include <array>
#include <iostream>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \brief The Shomate method to compute enthalpy and heat capacity
 * \tparam Scalar the type of the scalar values
 * \tparam intervals the static number of intervals for the temperature range, each interval has its own set of coefficients.
 *         If -1 the number of intervals is dynamic.
 */
template<class Scalar, int intervals = -1>
class ShomateMethod
{
    static_assert(intervals == -1 || intervals > 0, "Number of intervals must be -1 (dynamic) or greater than 0 (static).");

public:
    using CoefficientSet = struct { Scalar A, B, C, D, E, F, G, H; };
    using Temperatures = std::conditional_t<intervals == -1, std::vector<Scalar>, std::array<Scalar, std::size_t(intervals+1)>>;
    using Coefficients = std::conditional_t<intervals == -1, std::vector<CoefficientSet>, std::array<CoefficientSet, std::size_t(intervals)>>;

    /*!
     * \brief Constructor for Shomate method class
     * \param temperatures lower bound of the temperature intervals plus the upper bound of the last interval
     * \param coeffs contains the sets of coefficients for each temperature interval bounded by the entries of temperatures
     */
    constexpr ShomateMethod(const Temperatures& temperatures, const Coefficients& coeffs)
    : temperatures_(temperatures)
    , coeffs_(coeffs)
    {
        checkInput_();
    }

    /*!
     * \brief Return enthalpy in kJ/mol
     * \param temperature in K
     * \return Scalar
     */
    Scalar enthalpy(const Scalar temperature) const
    {
        const auto& p = paramsAtTemperature_(temperature);
        const Scalar t = temperature/1000.0;
        const Scalar standardEnthalpyDiff = p.A*t + p.B*t*t/2.0 + p.C*t*t*t/3.0 + p.D*t*t*t*t/4.0 - p.E/t + p.F - p.H;
        return standardEnthalpyDiff;
    }

    /*!
     * \brief Return heat capacity in J/(mol*K)
     * \param temperature in K
     * \return Scalar
     */
    Scalar heatCapacity(const Scalar temperature) const
    {
        const auto& p = paramsAtTemperature_(temperature);
        const Scalar t = temperature/1000.0;
        const Scalar heatCapacity = p.A + p.B*t + p.C*t*t + p.D*t*t*t + p.E/(t*t);
        return heatCapacity;
    }

private:
    const CoefficientSet& paramsAtTemperature_(const Scalar T) const
    {
        // check if T is smaller or higher than allowed min/max T
        if (T < temperatures_.front() || T > temperatures_.back())
        {
            if (!warningPrinted_)
            {
                std::cout << "Temperature "<< T << " [K] is out of range. Enthalpy values are extrapolated." << std::endl;
                warningPrinted_ = true;
            }
        }

        // find the interval for the given temperature
        const auto index = std::min<std::size_t>(
            coeffs_.size()-1,
            std::distance(
                temperatures_.begin(),
                std::lower_bound(temperatures_.begin(), temperatures_.end(), T)
            )
        );

        return coeffs_[index];
    }

    void checkInput_() const
    {
        if constexpr (intervals == -1)
        {
            if (temperatures_.size() < 2)
                DUNE_THROW(Dune::InvalidStateException, "Temperature range must have at least two entries.");

            for (size_t i = 0; i < temperatures_.size()-1; i++)
                if (temperatures_[i] >= temperatures_[i+1])
                    DUNE_THROW(Dune::InvalidStateException, "Temperature range must be strictly increasing.");

            if (temperatures_.size()-1 != coeffs_.size())
                DUNE_THROW(Dune::InvalidStateException, "If temperature vector is size n+1, there must be n coefficient sets.");
        }
    }

    Temperatures temperatures_;
    Coefficients coeffs_;
    static bool warningPrinted_;
};

template<class Scalar, int intervals>
bool ShomateMethod<Scalar, intervals>::warningPrinted_ = false;

} // namespace Dumux

#endif
