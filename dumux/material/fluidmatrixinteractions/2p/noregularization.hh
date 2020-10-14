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
 * \ingroup Fluidmatrixinteractions
 * \brief   A tag to turn off regularization and it's overhead
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_NO_REGULARIZATION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_NO_REGULARIZATION_HH

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A tag to turn off regularization and it's overhead
 */
struct NoRegularization
{
    //! Empty parameter structure
    template<class S> struct Params {};

    //! We are always equal to other instances of our kind
    bool operator== (const NoRegularization& o) const
    { return true; }
};

} // end namespace Dumux::FluidMatrix

#endif
