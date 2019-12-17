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
 * \brief Defines the properties required for (immiscible) two-phase sequential models.
 */
#ifndef DUMUX_2PPROPERTIES_SEQUENTIAL_HH
#define DUMUX_2PPROPERTIES_SEQUENTIAL_HH

// Dumux-includes
#include <dumux/porousmediumflow/sequential/properties.hh>
#include "indices.hh"
#include <dumux/material/spatialparams/sequentialfv.hh>

namespace Dumux {
////////////////////////////////
// Forward declarations
////////////////////////////////
template <class TypeTag, bool enableCompressibility>
class CellData2P;

////////////////////////////////
// Properties
////////////////////////////////
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The TypeTag for sequential two-phase problems
// Create new type tags
namespace TTag {
struct SequentialTwoP { using InheritsFrom = std::tuple<SequentialModel>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct  SaturationFormulation { using type = UndefinedProperty; }; //!< The formulation of the saturation model
template<class TypeTag, class MyTypeTag>
struct  VelocityFormulation { using type = UndefinedProperty; }; //!< The type of velocity reconstructed for the transport model
template<class TypeTag, class MyTypeTag>
struct  EnableCompressibility { using type = UndefinedProperty; };//!< Returns whether compressibility is allowed
} // end namespace Properties
} // end namespace Dumux

#include <dumux/porousmediumflow/sequential/variableclass.hh>
#include <dumux/porousmediumflow/2p/sequential/celldata.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidstates/isothermalimmiscible.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
//! Set number of equations to 2 for isothermal two-phase models
template<class TypeTag>
struct NumEq<TypeTag, TTag::SequentialTwoP> { static constexpr int value = 2; };

//! Set number of phases to 2 for two-phase models
template<class TypeTag>
struct NumPhases<TypeTag, TTag::SequentialTwoP> { static constexpr int value = 2; };//!< The number of phases in the 2p model is 2

//! Set number of components to 1 for immiscible two-phase models
template<class TypeTag>
struct NumComponents<TypeTag, TTag::SequentialTwoP> { static constexpr int value = 1; }; //!< Each phase consists of 1 pure component

//! Set \f$p_w\f$-\f$S_w\f$ formulation as default two-phase formulation
template<class TypeTag>
struct Formulation<TypeTag, TTag::SequentialTwoP> { static constexpr int value = SequentialTwoPCommonIndices::pwsw; };

//! Chose the set of indices depending on the chosen formulation
template<class TypeTag>
struct Indices<TypeTag, TTag::SequentialTwoP>
{
    using type = SequentialTwoPIndices<getPropValue<TypeTag, Properties::Formulation>(), 0>;
};

//! Set the default pressure formulation according to the chosen two-phase formulation
template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::SequentialTwoP> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::pressureType; };

//! Set the default saturation formulation according to the chosen two-phase formulation
template<class TypeTag>
struct SaturationFormulation<TypeTag, TTag::SequentialTwoP> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::saturationType; };

//! Set the default velocity formulation according to the chosen two-phase formulation
template<class TypeTag>
struct VelocityFormulation<TypeTag, TTag::SequentialTwoP> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::velocityDefault; };

//! Disable compressibility by default
template<class TypeTag>
struct EnableCompressibility<TypeTag, TTag::SequentialTwoP> { static constexpr bool value = false; };

//! Set general sequential VariableClass as default
template<class TypeTag>
struct Variables<TypeTag, TTag::SequentialTwoP> { using type = VariableClass<TypeTag>; };

//! Set standart CellData of immiscible two-phase models as default
template<class TypeTag>
struct CellData<TypeTag, TTag::SequentialTwoP> { using type = CellData2P<TypeTag, getPropValue<TypeTag, Properties::EnableCompressibility>()>; };

//! Set default fluid state
template<class TypeTag>
struct FluidState<TypeTag, TTag::SequentialTwoP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = IsothermalImmiscibleFluidState<Scalar, FluidSystem>;
};

//! The spatial parameters to be employed. Use SequentialFVSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SequentialTwoP> { using type = SequentialFVSpatialParams<TypeTag>; };
// \}
} // end namespace Properties
} // end namespace Dumux

#endif
