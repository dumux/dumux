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
 * \ingroup FreeflowModels
 * \brief This file contains different functions for estimating turbulence properties.
 */

#ifndef DUMUX_TURBULENCE_PROPERTIES_HH
#define DUMUX_TURBULENCE_PROPERTIES_HH

#include <iostream>

#include<dune/common/fvector.hh>

namespace Dumux {

/*!
 * \brief This class contains different functions for estimating turbulence properties.
 */
template<class Scalar, unsigned dim, bool verbose = false>
class TurbulenceProperties
{
public:
    /*!
     * \brief Estimates dimensionless wall distance \f$ y^+ \f$ based on a formula given in
     *        http://www.cfd-online.com/Wiki/Y_plus_wall_distance_estimation
     */
    Scalar yPlusEstimation(const Scalar velocity,
                           const Dune::FieldVector<Scalar, dim> position,
                           const Scalar kinematicViscosity,
                           const Scalar density,
                           int yCoordDim=dim-1,
                           bool print=verbose) const
    {
        using std::pow;
        using std::log10;
        using std::sqrt;
        const Scalar re_x = reynoldsNumber(velocity, position[0], kinematicViscosity, false);
        const Scalar c_f = pow((2.0 * log10(re_x) - 0.65), -2.3); // for re_x < 10^9
        const Scalar wallShearStress = 0.5 * c_f * density * velocity * velocity;
        const Scalar frictionVelocity = sqrt(wallShearStress / density);
        const Scalar yPlus = position[yCoordDim] * frictionVelocity / kinematicViscosity;
        if (print)
        {
            std::cout << "turbulence properties at (";
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                std::cout << position[dimIdx] << ",";
            std::cout << ")" << std::endl;
            std::cout << "estimated Re_x  : " << re_x << " [-]" << std::endl;
            std::cout << "estimated c_f   : " << c_f << " [-]" << std::endl;
            std::cout << "estimated tau_w : " << wallShearStress << " [kg/(m*s^2)]" << std::endl;
            std::cout << "estimated UStar : " << frictionVelocity << " [m/s]" << std::endl;
            std::cout << "estimated yPlus : " << yPlus << " [-]" << std::endl;
            std::cout << std::endl;
        }
        return yPlus;
    }

    /*!
     * \brief Estimates the entrance length for this pipe
     */
    Scalar entranceLength(const Scalar velocity,
                          const Scalar diameter,
                          const Scalar kinematicViscosity,
                          bool print=verbose) const
    {
        using std::pow;
        const Scalar re_d = reynoldsNumber(velocity, diameter, kinematicViscosity, false);
        const Scalar entranceLength = 4.4 * pow(re_d, 1.0/6.0) * diameter;
        if (print)
        {
            std::cout << "estimated Re_d  : " << re_d << " [-]" << std::endl;
            std::cout << "estimated l_ent : " << entranceLength << " [m]"<< std::endl;
            std::cout << std::endl;
        }
        return entranceLength;
    }

    /*!
     * \brief Calculates the Reynolds number
     */
    Scalar reynoldsNumber(const Scalar velocity,
                          const Scalar charLengthScale/*e.g. diameter*/,
                          const Scalar kinematicViscosity,
                          bool print=verbose) const
    {
        using std::abs;
        return abs(velocity * charLengthScale / kinematicViscosity);
    }

    /*!
     * \brief Estimates the turbulence intensity based on a formula given
     *        in the ANSYS Fluent user guide \cite ANSYSUserGuide12
     */
    Scalar turbulenceIntensity(const Scalar reynoldsNumber,
                               bool print=verbose) const
    {
        using std::pow;
        const Scalar turbulenceIntensity = 0.16 * pow(reynoldsNumber, -0.125);
        if (print)
        {
            std::cout << "estimated I     : " << turbulenceIntensity << " [-]" << std::endl;
        }
        return turbulenceIntensity;
    }

    /*!
     * \brief Estimates the turbulence length scale based on a formula given
     *        in the ANSYS Fluent user guide \cite ANSYSUserGuide12
     */
    Scalar turbulenceLengthScale(const Scalar charLengthScale/*e.g. diameter*/,
                                 bool print=verbose) const
    {
        const Scalar turbulenceLengthScale = 0.07 * charLengthScale;
        if (print)
        {
            std::cout << "estimated l_turb: " << turbulenceLengthScale << " [m]" << std::endl;
        }
        return turbulenceLengthScale;
    }

    /*!
     * \brief Estimates the turbulent kinetic energy based on a formula given
     *        in the ANSYS Fluent user guide \cite ANSYSUserGuide12
     */
    Scalar turbulentKineticEnergy(const Scalar velocity,
                                  const Scalar diameter,
                                  const Scalar kinematicViscosity,
                                  bool print=verbose) const
    {
        const Scalar re_d = reynoldsNumber(velocity, diameter, kinematicViscosity, false);
        const Scalar k = 1.5 * velocity * velocity
                         * turbulenceIntensity(re_d, false) * turbulenceIntensity(re_d, false);
        if (print)
        {
            std::cout << "estimated k     : " << k << " [m^2/s^2]" << std::endl;
        }
        return k;
    }

    /*!
     * \brief Estimates the dissipation based on a formula given
     *        in the ANSYS Fluent user guide \cite ANSYSUserGuide12
     */
    Scalar dissipation(const Scalar velocity,
                       const Scalar diameter,
                       const Scalar kinematicViscosity,
                       bool print=verbose) const
    {
        using std::pow;
        const Scalar k = turbulentKineticEnergy(velocity, diameter, kinematicViscosity, false);
        const Scalar factor = 0.1643; // = cMu^(3/4) = 0.09^(3/4)
        const Scalar epsilon = factor * pow(k, 1.5) / turbulenceLengthScale(diameter, false);
        if (print)
        {
            std::cout << "estimated eps.  : " << epsilon << " [m^2/s^3]" << std::endl;
        }
        return epsilon;
    }

    /*!
     * \brief Estimates the dissipation rate based on a formula given
     *        in the ANSYS Fluent user guide \cite ANSYSUserGuide12
     * \f[ \omega = \frac{k^{1/2}}{C_{\mu}^{1/4}L} \f]
     */
    Scalar dissipationRate(const Scalar velocity,
                           const Scalar diameter,
                           const Scalar kinematicViscosity,
                           bool print=verbose) const
    {
        using std::pow;
         const Scalar k = turbulentKineticEnergy(velocity, diameter, kinematicViscosity, false);
         const Scalar factor = 0.54772; // = cMu^(1/4) = 0.09^(1/4)
         const Scalar L = turbulenceLengthScale(diameter , false);
         const Scalar omega = pow(k, 0.5) / (L*factor);
        if (print)
        {
            std::cout << "estimated omega : " << omega << " [1/s]" << std::endl;

        }
        return omega;
    }

    /*!
     * \brief Estimates the viscosity tilde based on a formula given in
     *        in the ANSYS Fluent user guide \cite ANSYSUserGuide12
     */
    Scalar viscosityTilde(const Scalar velocity,
                          const Scalar diameter,
                          const Scalar kinematicViscosity,
                          bool print=verbose) const
    {
        using std::abs;
        using std::sqrt;
        const Scalar re_d = reynoldsNumber(velocity, diameter, kinematicViscosity, false);
        const Scalar viscosityTilde = sqrt(1.5) * abs(velocity)
                                      * turbulenceIntensity(re_d)
                                      * turbulenceLengthScale(diameter);
        if (print)
        {
            std::cout << "estimated nu~   : " << viscosityTilde << " [m^2/s]" << std::endl;
        }
        return viscosityTilde;
    }
};
} // end namespace Dumux

#endif
