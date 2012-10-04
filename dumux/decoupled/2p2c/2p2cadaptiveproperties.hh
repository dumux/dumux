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
 * \ingroup IMPEC
 * \ingroup IMPETProperties
 *
 * \file
 *
 * \brief Defines the properties required for the decoupled 2p2c models.
 */
#ifndef DUMUX_2P2CADAPTIVE_PROPERTIES_HH
#define DUMUX_2P2CADAPTIVE_PROPERTIES_HH

#include <dumux/decoupled/2p2c/2p2cproperties.hh>

namespace Dumux
{

//****** forward declarations  ******//
template<class TypeTag>
struct DecoupledTwoPTwoCIndicesAdaptive;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
//! The type tag for the compositional two-phase problems
NEW_TYPE_TAG(DecoupledTwoPTwoCAdaptive, INHERITS_FROM(DecoupledTwoPTwoC));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( GridAdaptEnableMultiPointFluxApproximation); //!< HangingNode: Two-point flux approximation (false) or mpfa (true)
NEW_PROP_TAG( GridAdaptEnableSecondHalfEdge ); //!< Uses second interaction volume for second half-edge in 2D
}}

//Dumux includes
#include <dumux/decoupled/2p2c/fvpressure2p2cadaptive.hh>
#include <dumux/decoupled/2p2c/fvtransport2p2cadaptive.hh>
#include <dumux/decoupled/2p2c/variableclass2p2cadaptive.hh>
#include <dumux/decoupled/2p2c/cellData2p2cadaptive.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
SET_BOOL_PROP(DecoupledTwoPTwoCAdaptive, AdaptiveGrid, true);
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, GridTypeIndices, GridTypes); //! Property not used but default necessary for mpfa2p
SET_BOOL_PROP(DecoupledTwoPTwoCAdaptive, GridAdaptEnableSecondHalfEdge, true); //!< Uses second interaction volume for second half-edge in 2D
SET_BOOL_PROP(DecoupledTwoPTwoCAdaptive, GridAdaptEnableMultiPointFluxApproximation, true); //!< applies an mpfa method around hanging nodes
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, CellData, CellData2P2CAdaptive<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, Variables, VariableClass2P2CAdaptive<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, Indices, DecoupledTwoPTwoCIndicesAdaptive<TypeTag>);
// Set the model properties
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, TransportModel, FVTransport2P2CAdaptive<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, PressureModel, FVPressure2P2CAdaptive<TypeTag>);
}


/*!
 * \brief Missing indices to the mpfa2p model.
 *
 * Compositional adaptive models use the 2p implementation with mpfa to
 * calculate the transmissibility (and nothing else). To couple both modules,
 * several Indice have to be present. Here, we apply dummy values to avoid
 * errors in case those Indice are really applied somewhere.
 */
template <class TypeTag>
struct DecoupledTwoPTwoCIndicesAdaptive : DecoupledTwoPTwoCIndices<TypeTag>
{
    static const int pressureIdx = 0;
    static const int saturationIdx = 0;
    static const int pressEqIdx = 0;
    static const int satEqIdx = 0;
};

// \}

}

#endif
