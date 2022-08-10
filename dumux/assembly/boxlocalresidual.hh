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
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief Calculates the element-wise residual for the box scheme
 */
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#warning "This header is deprecated and will be removed after 3.6"
#include "cvfelocalresidual.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup BoxDiscretization
 * \brief The element-wise residual for the box scheme
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
using BoxLocalResidual [[deprecated("Will be removed after 3.6")]] = CVFELocalResidual<TypeTag>;

} // end namespace Dumux

#endif
