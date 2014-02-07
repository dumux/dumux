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
 * \brief Defines the indices required for the compositional Stokes box model.
 */
#ifndef DUMUX_STOKES2C_INDICES_HH
#define DUMUX_STOKES2C_INDICES_HH

#include <dumux/freeflow/stokes/stokesindices.hh>

namespace Dumux
{
// \{

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the compositional Stokes box model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct Stokes2cCommonIndices : public StokesCommonIndices<TypeTag>
{
    // Phase index
    //! Index of the employed phase in case of a two-phase fluidsystem (set by default to nPhase)
    static const int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

    // Component indices
    static const int phaseCompIdx = phaseIdx; //!< The index of the main component of the considered phase
    //! The index of the transported (minor) component; ASSUMES phase indices of 0 and 1
    static const int transportCompIdx = (unsigned int)(1-phaseIdx);

    // equation and primary variable indices
    static const int dim = StokesCommonIndices<TypeTag>::dim;
    static const int transportEqIdx = PVOffset + dim+1; //!< The index for the transport equation

    //! The index of the mass or mole fraction of the transported component in primary variable vectors
    static const int massOrMoleFracIdx = transportEqIdx;
};
} // end namespace

#endif
