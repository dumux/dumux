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
 * \ingroup ThreePThreeCTests
 * \brief The properties of the non-isothermal gas injection problem where a gas (e.g. steam/air) is
 *        injected into a unsaturated porous medium with a residually trapped NAPL contamination.
 */

#ifndef DUMUX_KUEVETTE3P3CNIPROPERTIES_HH
#define DUMUX_KUEVETTE3P3CNIPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p3c/model.hh>

#define ISOTHERMAL 0

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct Kuevette { using InheritsFrom = std::tuple<ThreePThreeCNI>; };
struct KuevetteBox { using InheritsFrom = std::tuple<Kuevette, BoxModel>; };
struct KuevetteCCTpfa { using InheritsFrom = std::tuple<Kuevette, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Kuevette> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Kuevette> { using type = KuevetteProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Kuevette>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = KuevetteSpatialParams<GridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Kuevette>
{ using type = FluidSystems::H2OAirMesitylene<GetPropType<TypeTag, Properties::Scalar>>; };
} // end namespace Dumux::Properties

#endif
