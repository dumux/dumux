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
 *
 * \brief Specification of the material parameters
 *       for the PNM constitutive relations.
 */
#ifndef DUMUX_PNM_LOCAL_RULES_PARAMS_HH
#define DUMUX_PNM_LOCAL_RULES_PARAMS_HH

namespace Dumux
{

/*!
 * \brief Specification of the material parameters
 *       for the PNM constitutive relations.
 *
 *        \ingroup fluidmatrixinteractionsparams
 *
 */
template <class Scalar>
class PNMLocalRulesParams
{
public:

    PNMLocalRulesParams(const Scalar surfaceTension,
                        const Scalar contactAngle,
                        const Scalar inscribedRadius)
    : surfaceTension_(surfaceTension)
    , contactAngle_(contactAngle)
    , inscribedRadius_(inscribedRadius)
    {}

    /*!
    * \brief Returns the inscribed pore radius \f$\mathrm{[m]}\f$
    */
    Scalar poreRadius() const noexcept
    { return inscribedRadius_;}

    /*!
    * \brief Returns the contact angle \f$\mathrm{[rad]}\f$
    */
    Scalar contactAngle() const noexcept
    { return contactAngle_;}

     /*!
     * \brief Returns the interfacial tension \f$\mathrm{[kg/s^2]}\f$
     */
    Scalar surfaceTension() const noexcept
    { return surfaceTension_; }

private:

    Scalar surfaceTension_ ;
    Scalar contactAngle_ ;
    Scalar inscribedRadius_;

};
} // namespace Dumux

#endif
