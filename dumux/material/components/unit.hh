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
 * \ingroup Components
 * \brief A component using a value of 1 for all fluid properties.
 */
#ifndef DUMUX_UNIT_HH
#define DUMUX_UNIT_HH

#include <dune/common/deprecated.hh>
#include "component.hh"

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief A component using a value of one for all fluid properties.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class DUNE_DEPRECATED_MSG("Use Components::Constant<id, Scalar> instead. The default is a unit fluid system.")
Unit : public Component<Scalar, Unit<Scalar> >
{

public:
    /*!
     * \brief A human readable name for the component.
     */
    static std::string name()
    { return "Unit"; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Sets the density to 1 \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1.0;
    }

    /*!
     * \brief Sets the viscosity to 1 \f$\mathrm{[Pa*s]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1.0;
    }

};

} // end namespace

#endif
