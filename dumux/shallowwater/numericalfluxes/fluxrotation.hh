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
 * \ingroup SweModel
 * \brief Rotate the flux using the rotation invarianz of the SWEs
 *        two-point flux approximation.
 */
#ifndef DUMUX_SHALLOWWATER_NUMERICALFLUXES_FLUXROTATION_HH
#define DUMUX_SHALLOWWATER_NUMERICALFLUXES_FLUXROTATION_HH


#include <dumux/common/math.hh>

namespace Dumux
{

inline void stateRotation(const auto& nxy, auto& cellStates)
  {
    //Rotate variables u.v before hllc call
    auto temp = cellStates[1];
    cellStates[1] =  nxy[0] * temp + nxy[1] * cellStates[2];
    cellStates[2] = -nxy[1] * temp + nxy[0] * cellStates[2];
  }


inline void rotateFluxBack(const auto& nxy, auto& flux)
  {
    //rotate the fluxes back after hllc call
    auto temp = flux[1];
    flux[1] = nxy[0] * temp - nxy[1] * flux[2];
    flux[2] = nxy[1] * temp + nxy[0] * flux[2];
  }


} // end namespace Dumux

#endif
