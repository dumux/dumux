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
 * \ingroup CCMpfaDiscretization
 * \brief Stores the face indices corresponding to the neighbors of an element
 *        that contribute to the derivative calculation. Depending on if an
 *        mpfa scheme leads to a symmetric/unsymmetric sparsity pattern, the
 *        adequate implementation of the connectiviy map is chosen.
 */
#ifndef DUMUX_CC_MPFA_CONNECTIVITY_MAP_HH
#define DUMUX_CC_MPFA_CONNECTIVITY_MAP_HH

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/connectivitymap.hh>

namespace Dumux {

//! Forward declaration of method specific implementation of the assembly map
template<class GridGeometry, MpfaMethods method>
class CCMpfaConnectivityMap;

//! The o-method can use the simple (symmetric) assembly map
template<class GridGeometry>
class CCMpfaConnectivityMap<GridGeometry, MpfaMethods::oMethod> : public CCSimpleConnectivityMap<GridGeometry> {};
} // end namespace Dumux

#endif
