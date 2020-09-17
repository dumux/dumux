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
 * \brief Implementation of the capillary pressure and
 * relative permeability <-> saturation relations according to Joekar-Niasar et al., 2010.
 *
 */
#ifndef DUMUX_PNM_2P_BASE_LOCAL_RULES_HH
#define DUMUX_PNM_2P_BASE_LOCAL_RULES_HH

#include <dumux/porenetworkflow/common/poreproperties.hh>

namespace Dumux
{

struct TwoPLocalRulesBase
{
    /*!
     * \brief The parameter type
     * \tparam Scalar The scalar type
     */
    template<class Scalar>
    struct Params
    {
        Scalar poreRadius, contactAngle, surfaceTension;
        Pore::Shape shape;
    };

    template<class... Args>
    static double krw(Args&&...)
    { return 1.0; }

    template<class... Args>
    static double krn(Args&&...)
    { return 1.0; }

    template<class... Args>
    static double dkrw_dsw(Args&&...)
    { return 0.0; }

    template<class... Args>
    static double dkrn_dsw(Args&&...)
    { return 0.0; }

};

}

#endif // DUMUX_PNM_LOCAL_RULES_HH
