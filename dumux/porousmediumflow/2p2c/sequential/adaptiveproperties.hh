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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup SequentialTwoPTwoCModel
 * \brief Defines the properties required for the adaptive sequential 2p2c models.
 */
#ifndef DUMUX_2P2CADAPTIVE_PROPERTIES_HH
#define DUMUX_2P2CADAPTIVE_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p2c/sequential/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>

namespace Dumux {

//****** forward declarations  ******//
template<class TypeTag>
struct SequentialTwoPTwoCIndicesAdaptive;
template<class TypeTag>
class FV2dPressure2P2CAdaptive;
template<class TypeTag>
class FV2dTransport2P2CAdaptive;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
//! The type tag for the compositional two-phase problems
// Create new type tags
namespace TTag {
struct SequentialTwoPTwoCAdaptive { using InheritsFrom = std::tuple<SequentialTwoPTwoC>; };
} // end namespace TTag
}}

//Dumux includes
#include <dumux/porousmediumflow/2p2c/sequential/fv2dpressureadaptive.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fv3dpressureadaptive.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fv2dtransportadaptive.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fv3dtransportadaptive.hh>
#include <dumux/porousmediumflow/2p2c/sequential/variableclassadaptive.hh>
#include <dumux/porousmediumflow/2p2c/sequential/celldataadaptive.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
template<class TypeTag>
struct AdaptiveGrid<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { static constexpr bool value = true; };
template<class TypeTag>
struct GridTypeIndices<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = GridTypes; }; //!< Property not used but default necessary for mpfa2p

template<class TypeTag>
struct CellData<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = CellData2P2CAdaptive<TypeTag>; };
template<class TypeTag>
struct Variables<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = VariableClass2P2CAdaptive<TypeTag>; };
template<class TypeTag>
struct Indices<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = SequentialTwoPTwoCIndicesAdaptive<TypeTag>; };
// Set the model properties
template<class TypeTag>
struct TransportModel<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = FV2dTransport2P2CAdaptive<TypeTag>; };
template<class TypeTag>
struct PressureModel<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = FV2dPressure2P2CAdaptive<TypeTag>; };
}


/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Missing indices to the mpfa2p model.
 *
 * Compositional adaptive models use the 2p implementation with mpfa to
 * calculate the transmissibility (and nothing else). To couple both modules,
 * several Indice have to be present. Here, we apply dummy values to avoid
 * errors in case those Indice are really applied somewhere.
 */
template <class TypeTag>
struct SequentialTwoPTwoCIndicesAdaptive : public SequentialTwoPTwoCIndices<TypeTag>
{
    static const int pressureIdx = 0;
    static const int saturationIdx = 0;
    static const int pressureEqIdx = 0;
    static const int satEqIdx = 0;
};

// \}

} // end namespace Dumux

#endif
