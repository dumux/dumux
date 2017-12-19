// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Binary coefficients for water and a ficticious component implemented in tutorial exercise 3a.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_MYCOMPRESSIBLECOMPONENT_HH
#define DUMUX_BINARY_COEFF_H2O_MYCOMPRESSIBLECOMPONENT_HH

namespace Dumux
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and a ficticious component implemented in tutorial exercise 3a
 *        The implementation of the missing methods in this file is part of exercise 3b.
 */
class H2O_MyCompressibleComponent
{
public:
    /*!
     * \brief Henry coefficient \f$[N/m^2]\f$ for the fictitous component in liquid water.
     */
    template <class Scalar>
    static Scalar henryMyCompressibleComponentInWater(Scalar temperature)
    {
        Scalar dumuxH = 1.5e-1 / 101.325; // unit [(mol/m^3)/Pa]
        dumuxH *= 18.02e-6;  //multiplied by molar volume of reference phase = water
        return 1.0/dumuxH; // [Pa]
    }

    /*!
     * \brief Henry coefficient \f$[N/m^2]\f$ for water in the ficticious component.
     */
    template <class Scalar>
    static Scalar henryWaterInMyCompressibleComponent(Scalar temperature)
    {
        // arbitrary
        return 1.0e8; // [Pa]
    }

    /*!
     * \brief Diffusion coefficient [m^2/s] for my ficticious component in liquid water or vice versa.
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        // arbitrary
        return 1.e-9;
    }
};

}
} // end namespace

#endif
