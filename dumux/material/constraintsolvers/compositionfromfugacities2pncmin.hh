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
 *
 * \brief Determines the fluid composition given the component
 *        fugacities and an arbitary equation of state.
 */
#ifndef DUMUX_COMPOSITION_FROM_FUGACITIES_2PNCMIN_HH
#define DUMUX_COMPOSITION_FROM_FUGACITIES_2PNCMIN_HH

#include "compositionfromfugacities.hh"
#include <dune/common/deprecated.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {
/*!
 * \ingroup ConstraintSolver
 * \brief Calculates the chemical equilibrium from the component
 *        fugacities in a phase.
 */
template <class Scalar, class FluidSystem>
class DUNE_DEPRECATED_MSG("CompositionFromFugacities2pncmin is deprecated. Use CompositionFromFugacities instead.")
  CompositionFromFugacities2pncmin
  : public CompositionFromFugacities<Scalar, FluidSystem>
{ };

template <class Scalar, class FluidSystem>
class DUNE_DEPRECATED_MSG("compositionFromFugacities2pncmin is deprecated. Use CompositionFromFugacities2pncmin (capital C) instead.")
  compositionFromFugacities2pncmin
  : public CompositionFromFugacities2pncmin<Scalar, FluidSystem>
{ };
} // end namespace Dumux

#endif
