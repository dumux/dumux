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
 * \ingroup Common
 * \file
 * \brief Defines capabilities recognized for problems in Dumux
 * \note Specialize the subsequent capabilities for your problem to optimize
 *       your program's effeciency.
 */
#ifndef DUMUX_CAPABILITIES_HH
#define DUMUX_CAPABILITIES_HH

#include <dune/common/deprecated.hh>

namespace Dumux
{

namespace Capabilities
{

//! If a problem is stationary (not time-dependent)
template <class Problem>
struct DUNE_DEPRECATED_MSG("isStationary is deprecated and will be removed!") isStationary
{
    //! by default all problems are instationary
    static const bool value = false;
};

} // namespace Capabilities

} // namespace Dumux

#endif
