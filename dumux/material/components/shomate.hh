// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Apply Shomate equations for enthalpy and heat capacity.
 *
 * The Shomate equations were inroduced in
 * "A method for evaluating and correlating thermodynamic data" by C. Howard Shomate in 1954, \cite Shomate1954.
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
#include <vector>
#include <stdexcept>
#include <iostream>

namespace Dumux {
template<class Scalar>
class ShomateMethod
{
public:
    /*!
     *  \brief Constructor for Shomate method class
     *
     * \param tempRange contains the temperature boundaries of the intervals defined by the Shomate method, the number of intervals can vary
     * \param coeffs contains the sets of coefficients for each temperature interval bounded by the entries of tempRange
     *
     */
    ShomateMethod(std::vector<Scalar> tempRange, std::vector<std::array<Scalar,8>> coeffs)
    : tempRange_(tempRange)
    , coeffs_(coeffs)
    {
        if (tempRange_.size()-1 != coeffs_.size())
        {
            DUNE_THROW(Dune::InvalidStateException,"Temperature range and coefficients must have the same size.");
        }
        for(size_t i=0; i<coeffs_.size(); i++)
        {
            if (coeffs_[i].size() != 8)
            {
                DUNE_THROW(Dune::InvalidStateException,"Shomate coeffs must have 8 elements.");
            }
        }
    }

    std::array<Scalar,8> paramsAtTemperature(const Scalar& T) const
    {
        // check if T is smaller or higher than allowed min/max T
        if (T < tempRange_[0])
        {
            if (!warningPrinted_)
            {
                std::cerr << "Temperature "<< T << " [K] is out of range. Enthalpy values are extrapolated." << std::endl;
                warningPrinted_ = true;
            }
            return coeffs_[0];
        }

        for (size_t i = 0; i < tempRange_.size()-1; i++)
        {
            if (T >= tempRange_[i] && T < tempRange_[i+1])
            {
                return coeffs_[i];
            }
        }

        if (!warningPrinted_)
        {
            std::cerr << "Temperature "<< T << " [K] is out of range. Enthalpy values are extrapolated." << std::endl;
            warningPrinted_ = true;
        }
        return coeffs_[tempRange_.size()-1];
    }

    /**
     * @brief return enthalpy in kJ/mol
     *
     * @param temperature in K
     * @param pressure in Pa
     * @return Scalar
     */
    Scalar enthalpy(const Scalar& temperature,
                       const Scalar& pressure) const
    {
        const auto& correctCoeffs = paramsAtTemperature(temperature);

        Scalar A = correctCoeffs[0];
        Scalar B = correctCoeffs[1];
        Scalar C = correctCoeffs[2];
        Scalar D = correctCoeffs[3];
        Scalar E = correctCoeffs[4];
        Scalar F = correctCoeffs[5];
        Scalar H = correctCoeffs[7];
        //convert temperature to correct format
        Scalar t = temperature/1000.0;
        //calculate standard enthalpy difference according to Shomate
        Scalar standardEnthalpyDiff = A*t + B*t*t/2.0 + C*t*t*t/3.0 + D*t*t*t*t/4.0 - E/t + F - H;
        //unit is kJ/mol
        return standardEnthalpyDiff;
    }

    /**
     * @brief return heat capacity in J/(mol*K)
     *
     * @param temperature in K
     * @param pressure in Pa
     * @return Scalar
     */
    Scalar heatCapacity(const Scalar& temperature,
                        const Scalar& pressure) const
    {
        auto correctCoeffs = paramsAtTemperature(temperature);
        Scalar A = correctCoeffs[0];
        Scalar B = correctCoeffs[1];
        Scalar C = correctCoeffs[2];
        Scalar D = correctCoeffs[3];
        Scalar E = correctCoeffs[4];
        Scalar t = temperature/1000.0;
        Scalar heatCapacity = A + B*t + C*t*t + D*t*t*t + E/(t*t);
        return heatCapacity;
    }

private:
    std::vector<Scalar> tempRange_;
    std::vector<std::array<Scalar,8>> coeffs_;
    static bool warningPrinted_;
};

template<class Scalar>
bool ShomateMethod<Scalar>::warningPrinted_ = false;
} // namespace Dumux
#endif // DUMUX_SHOMATE_HH
