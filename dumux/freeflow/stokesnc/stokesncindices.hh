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
 * \brief Defines the indices required for the compositional n component Stokes box model.
 */
#ifndef DUMUX_STOKESNC_INDICES_HH
#define DUMUX_STOKESNC_INDICES_HH

#include "stokesncproperties.hh"
#include <dumux/freeflow/stokes/stokesindices.hh>

namespace Dumux
{
// \{

/*!
 * \ingroup BoxStokesncModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the compositional n component Stokes box model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct StokesncCommonIndices : public StokesCommonIndices<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:

    // Dimension (copied for convenience)
    static const int dim = StokesCommonIndices<TypeTag>::dim; //!< Number of dimensions

    // Phase Index
    //! Index of the employed phase in case of a two-phase fluidsystem (set by default to nPhase)
    static const int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

    // Number of Components
    //! Number of components in employed fluidsystem
    static const int numComponents = FluidSystem::numComponents;

    // Component indices
    //! The index of the main component of the considered phase
    static const int phaseCompIdx = phaseIdx;
    //! The index of the first transported component; ASSUMES phase indices of 0 and 1
    static const int transportCompIdx = (unsigned int)(1-phaseIdx);

    // Transport equation indices
    //! The index of the mass conservation equation of the first component. In analogy to porous media models "conti"
    //! is used here to describe mass conservation equations, i.e total mass balance and transport equations.
    static const int conti0EqIdx = PVOffset + dim;
    //! The index of the mass balance equation sits on the slot of the employed phase
    static const int massBalanceIdx = conti0EqIdx + phaseCompIdx;
    //! The index of the transport equation for a two component model.
    //! For n>2 please don't use this index, because it looses its actual meaning.
    static const int transportEqIdx = conti0EqIdx + transportCompIdx;

    // Primary variables
    //! The index of the first mass or mole fraction of the transported component in primary variable vectors
    static const int massOrMoleFracIdx = transportEqIdx;
    //! The index of the pressure in primary variable vectors
    static const int pressureIdx = massBalanceIdx;
};
} // end namespace

#endif
