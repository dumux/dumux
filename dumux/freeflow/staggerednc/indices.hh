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
 * \brief  Defines the indices for the one-phase fully implicit model.
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_NC_INDICES_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_NC_INDICES_HH

#include <dumux/freeflow/staggered/indices.hh>

namespace Dumux
{
// \{
/*!
 * \ingroup NavierStokesModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the isothermal stokes model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct NavierStokesNCIndices : public NavierStokesCommonIndices<TypeTag, PVOffset>
{

};

// \}
} // end namespace

#endif
