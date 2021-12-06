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
 * \ingroup TwoPNCNCPTests
 * \brief The properties of the problem where liquid water is injected which has to go around an
 * TwoPNCSalt with \f$10^3\f$ lower permeability.
 */
#ifndef DUMUX_TEST_2PNC_NCP_PROPERTIES_HH
#define DUMUX_TEST_2PNC_NCP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/2pnc_ncp/model.hh>
#include <dumux/material/fluidsystems/brineair.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TwoPNCSalt { using InheritsFrom = std::tuple<TwoPNCNCP>; };
struct TwoPNCSaltBox { using InheritsFrom = std::tuple<TwoPNCSalt, BoxModel>; };
struct TwoPNCSaltCC { using InheritsFrom = std::tuple<TwoPNCSalt, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPNCSalt> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPNCSalt> { using type = TwoPNCSaltProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPNCSalt>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPNCSaltSpatialParams<GridGeometry, Scalar>;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPNCSalt>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BrineAir<Scalar, Components::H2O<Scalar>>;
};

// replace the main component balance eq with a total balance eq
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::TwoPNCSalt> { static constexpr int value = 10; };

} // end namespace Dumux::Properties

#endif
