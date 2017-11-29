// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *\ingroup Components
 * \brief  Rough estimate of some wetting phase for testing purposes.
 */
#ifndef DUMUX_SPE10WETTING_HH
#define DUMUX_SPE10WETTING_HH


#include <dumux/material/components/component.hh>

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief Rough estimate for testing purposes of some wetting phase.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SPE10Wetting : public Component<Scalar, SPE10Wetting<Scalar> >
{

public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "SPE10Wetting"; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Rough estimate of the density of oil \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 1025.18;
    }

    /*!
     * \brief Rough estimate of the viscosity of oil in \f$\mathrm{[Pa*s]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 0.0003;
    };

};

} // end namespace

#endif
