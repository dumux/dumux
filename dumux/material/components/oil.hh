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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *\ingroup Components
 * \brief  Rough estimate of some oil for testing purposes.
 */
#ifndef DUMUX_OIL_HH
#define DUMUX_OIL_HH

#warning use material/components/lnapl.hh instead

#include "lnapl.hh"

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief Rough estimate for testing purposes of some oil.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Oil : public LNAPL<Scalar>
{

public:
    DUNE_DEPRECATED_MSG("use LNAPL instead")
    Oil(){}
    /*!
     * \brief A human readable name for the water.
     */
    static const char *name()
    { return "Oil"; }
};

} // end namepace

#endif
