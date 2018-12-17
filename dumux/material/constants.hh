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
 * \ingroup Material
 * \brief A central place for various physical constants occuring in
 *        some equations.
 */
#ifndef DUMUX_CONSTANTS_HH
#define DUMUX_CONSTANTS_HH

#include <cmath>

namespace Dumux {

/*!
 * \ingroup Material
 * \brief A central place for various physical constants occuring in
 *        some equations.
 */
template<class Scalar>
class Constants
{
public:
    /*!
     * \brief The ideal gas constant \f$\mathrm{[J/(mol K)]}\f$
     */
    static constexpr Scalar R = 8.314472;

    /*!
     * \brief The Avogadro constant \f$\mathrm{[1/mol]}\f$
     */
    static constexpr Scalar Na = 6.02214179e23;

    /*!
     * \brief The Boltzmann constant \f$\mathrm{[J/K]}\f$
     */
    static constexpr Scalar kb = R / Na;

    /*!
     * \brief Speed of light in vacuum \f$\mathrm{[m/s]}\f$
     */
    static constexpr Scalar c = 299792458;

    /*!
     * \brief Faraday constant \f$\mathrm{[C/mol]}\f$
     *
     * Source: CODATA 2010
     */
    static constexpr Scalar F = 96485.3365;

    /*!
     * \brief Newtonian constant of gravitation \f$\mathrm{[m^3/(kg s^2)]}\f$
     */
    static constexpr Scalar G = 6.67428e-11;

    /*!
     * \brief Planck constant \f$\mathrm{[J s]}\f$
     */
    static constexpr Scalar h = 6.62606896e-34;

    /*!
     * \brief Reduced Planck constant \f$\mathrm{[J s]}\f$
     */
    static constexpr Scalar hRed = h / (2 * M_PI);
};

} // end namespace Dumux

#endif
