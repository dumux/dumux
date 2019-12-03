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
 * \ingroup SequentialTwoPModel
 * \brief Specifies the properties for immiscible 2p transport
 */
#ifndef DUMUX_FV_TRANSPORT_PROPERTIES_2P_HH
#define DUMUX_FV_TRANSPORT_PROPERTIES_2P_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/properties.hh>

namespace Dumux {
namespace Properties {
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for two-phase problems using a standard finite volume model
// Create new type tags
namespace TTag {
struct FVTransportTwoP { using InheritsFrom = std::tuple<TransportTwoP>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//! Bool property which tells the transport model if it should use constitutive relations which
//! are precomputed at the begin of the time step or if it should recompute the relations
template<class TypeTag, class MyTypeTag>
struct  PrecomputedConstRels  { using type = UndefinedProperty; };
} // end namespace Properties
} // end namespace Dumux

#include "evalcflfluxdefault.hh"
#include "dumux/porousmediumflow/2p/sequential/transport/cellcentered/diffusivepart.hh"
#include "dumux/porousmediumflow/2p/sequential/transport/cellcentered/convectivepart.hh"
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/saturation.hh>

namespace Dumux {
namespace Properties {
//! Set the default implementation of the cfl-condition
template<class TypeTag>
struct EvalCflFluxFunction<TypeTag, TTag::FVTransportTwoP> { using type = EvalCflFluxDefault<TypeTag>; };
//! Set the default implementation of a diffusive flux -> diffusive flux dissabled
template<class TypeTag>
struct CapillaryFlux<TypeTag, TTag::FVTransportTwoP> { using type = DiffusivePart<TypeTag>; };
//! Set the default implementation of an additional convective flux -> additional convective flux dissabled
template<class TypeTag>
struct GravityFlux<TypeTag, TTag::FVTransportTwoP> { using type = ConvectivePart<TypeTag>; };
//! Set PrecomputedConstRels flag <tt>true</tt> as default
template<class TypeTag>
struct PrecomputedConstRels<TypeTag, TTag:: FVTransportTwoP> { static constexpr bool value = true; };
//! Set finite volume implementation of the two-phase saturation equation as default saturation model
template<class TypeTag>
struct TransportModel<TypeTag, TTag::FVTransportTwoP> { using type = FVSaturation2P<TypeTag>; };
} // end namespace Properties
} // end namespace Dumux

#endif
