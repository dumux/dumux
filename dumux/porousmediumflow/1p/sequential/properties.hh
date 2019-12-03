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
 * \ingroup SequentialOnePModel
 * \brief Defines the properties required for the single phase sequential model.
 */

#ifndef DUMUX_1PPROPERTIES_HH
#define DUMUX_1PPROPERTIES_HH

//Dumux-includes
#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/material/spatialparams/sequentialfv1p.hh>

namespace Dumux {

////////////////////////////////
// forward declarations
////////////////////////////////
template <class TypeTag>
class CellData1P;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase problem
// Create new type tags
namespace TTag {
struct SequentialOneP { using InheritsFrom = std::tuple<SequentialModel>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct  Fluid  { using type = UndefinedProperty; };          // The fluid for one-phase models
} // end namespace Properties
} // end namespace Dumux

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/porousmediumflow/sequential/variableclass.hh>
#include <dumux/porousmediumflow/1p/sequential/celldata.hh>
#include <dumux/porousmediumflow/1p/sequential/indices.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Set number of equations to 1 for isothermal one-phase models
template<class TypeTag>
struct NumEq<TypeTag, TTag::SequentialOneP> { static constexpr int value = 1; };

//! Set number of phases to 1 for one-phase models
template<class TypeTag>
struct NumPhases<TypeTag, TTag::SequentialOneP> { static constexpr int value = 1; };

//! Each phase consists of 1 pure component
template<class TypeTag>
struct NumComponents<TypeTag, TTag::SequentialOneP> { static constexpr int value = 1; };

//! Chose the set of indices for the one-phase formulation
template<class TypeTag>
struct Indices<TypeTag, TTag::SequentialOneP> { using type = SequentialOnePCommonIndices; };

//! Set general sequential VariableClass as default
template<class TypeTag>
struct Variables<TypeTag, TTag::SequentialOneP> { using type = VariableClass<TypeTag>; };

//! Set standart CellData of immiscible one-phase models as default
template<class TypeTag>
struct CellData<TypeTag, TTag::SequentialOneP> { using type = CellData1P<TypeTag>; };

//! The spatial parameters to be employed. Use BoxSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SequentialOneP> { using type = SequentialFVSpatialParamsOneP<TypeTag>; };
} // end namespace Properties
} // end namespace Dumux
#endif
