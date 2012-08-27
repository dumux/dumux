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
 *
 * \brief Defines the indices required for the non-isothermal compositional Stokes box model.
 */
#ifndef DUMUX_STOKES2CNI_INDICES_HH
#define DUMUX_STOKES2CNI_INDICES_HH

#include <dumux/freeflow/stokes2c/stokes2cindices.hh>

namespace Dumux
{
// \{

/*!
 * \ingroup BoxStokes2cniModel
 * \ingroup BoxIndices
 * \brief Enumerations for the non-isothermal compositional Stokes model
 */
template <class TypeTag, int PVOffset=0>
struct Stokes2cniCommonIndices : public Stokes2cCommonIndices<TypeTag, PVOffset>
{
public:
    // number of dimensions
    static const int dim = StokesCommonIndices<TypeTag>::dim;
    static const int energyEqIdx = PVOffset + dim+2; //! The index for the energy balance equation
    static const int energyIdx = energyEqIdx; // \deprecated Use energyEqIdx instaed
    static const int temperatureIdx = energyEqIdx; //! The index for temperature in primary variable vectors
};
} // end namespace

#endif
