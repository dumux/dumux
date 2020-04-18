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

#ifndef DUMUX_EXAMPLE_TWOP_INFILTRATION_PROPERTIES_HH
#define DUMUX_EXAMPLE_TWOP_INFILTRATION_PROPERTIES_HH

// ## The file `properties.hh`
// [[content]]
//
// ### Includes
// The header includes will be mentioned in the text below.
// [[details]] header includes
#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/2p/model.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include "spatialparams.hh"
#include "problem.hh"
// [[/details]]
//
// ### Type tag definition
// All properties are defined in the (nested) namespace
// `Dumux::Properties`. To get and set properties, we need the definitions and implementations from the
// header `dumux/common/properties.hh` included above.
namespace Dumux::Properties {

// First, a so-called TypeTag is created. Properties are traits specialized for this TypeTag (a simple `struct`).
//
// >>>
// :white_check_mark: The properties of other TypeTags are inherited if the alias `InheritsFrom` is present.
// These other TypeTags are listed in form of a `std::tuple` in order of importance.
// >>>
//
// Here, properties from the two-phase flow model (`TTag::Twop`) and the
// cell-centered finite volume scheme with two-point-flux approximation (`TTag::CCTpfaModel`)
// are inherited. These other TypeTag definitions can be found in the included
// headers `dumux/porousmediumflow/2p/model.hh` and `dumux/discretization/cctpfa.hh`.
namespace TTag {
struct PointSourceExample { using InheritsFrom = std::tuple<TwoP, CCTpfaModel>; };
}
// ### Property specializations
// Next, we specialize the properties `Problem` and `SpatialParams` for our new TypeTag and
// set the type to our problem and spatial parameter classes implemented
// in `problem.hh` and `spatialparams.hh`.
// [[codeblock]]
template<class TypeTag>
struct Problem<TypeTag, TTag::PointSourceExample>
{ using type = PointSourceProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PointSourceExample>
{
    // two local aliases for convenience and readability
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using type = TwoPTestSpatialParams<GridGeometry, Scalar>;
};
// [[/codeblock]]

// The `Grid` property tells the
// simulator to use ALUGrid - an unstructured grid manager - here
// configured for grid and coordinate dimensions `2`,
// hexahedral element types (`Dune::cube`) and non-conforming refinement mode.
// `Dune::ALUGrid` is declared in the included header `dune/alugrid/grid.hh`
// from the Dune module `dune-alugrid`.
template<class TypeTag>
struct Grid<TypeTag, TTag::PointSourceExample>
{ using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };

// The `FluidSystem` property specifies which fluids are used.
// This fluid system is composed of two immiscible liquid phases which are made up
// entirely of its respective main components `SimpleH2O` (a water component with constant properties)
// and `Trichloroethene` (a DNAPL). The components, phases, and the fluid system are implemented in
// the headers `dumux/material/components/simpleh2o.hh`,
// `dumux/material/components/trichloroethene.hh`,
// `dumux/material/fluidsystems/1pliquid.hh`,
// `dumux/material/fluidsystems/2pimmiscible.hh`
// included above.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PointSourceExample>
{
  using Scalar = GetPropType<TypeTag, Properties::Scalar>;
  using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
  using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;

  using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// The two-phase model implements different primary variable formulations.
// Here we choose the pressure of the first phase and the saturation of the second phase.
// The order of phases is specified by the fluid system.
// In this case that means that the primary variables are water pressure and DNAPL saturation.
template<class TypeTag>
struct Formulation<TypeTag, TTag::PointSourceExample>
{ static constexpr auto value = TwoPFormulation::p0s1; };

} // end namespace Dumux::Properties
// [[/content]]
#endif
