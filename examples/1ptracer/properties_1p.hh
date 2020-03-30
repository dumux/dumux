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

#ifndef DUMUX_ONEP_TEST_PROPERTIES_HH
#define DUMUX_ONEP_TEST_PROPERTIES_HH

// ## Property definitions (`properties_1p.hh`)
//
// This file defines the `TypeTag` used for the single-phase flow simulation,
// for which we then define the necessary properties.
//
// ### Include files
// <details>
// The `TypeTag` defined for this simulation will inherit all properties from the
// `OneP` type tag, a convenience type tag that predefines most of the required
// properties for single-phase flow simulations in DuMuX.
#include <dumux/porousmediumflow/1p/model.hh>

// We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>
// In this example, we want to discretize the equations with the cell centered finite volume
// scheme using two-point-flux approximation:
#include <dumux/discretization/cctpfa.hh>
// The fluid properties are specified in the following headers (we use liquid water as the fluid phase):
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// The local residual for incompressible flow is included.
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

// We include the problem and spatial parameters headers used for this simulation.
#include "problem_1p.hh"
#include "spatialparams_1p.hh"
// </details>
//
// ### Definition of the `TypeTag` used for the single-phase problem
//
// A `TypeTag` for our simulation is created which inherits from the one-phase flow model
// (defined by the type tag `OneP`) and the cell centered finite volume scheme using two-point-flux
// approximation (defined by the type tag `CCTpfaModel`).
// The `OneP` type tag predefines most of the required properties for single-phase flow simulations
// in DuMuX. However, some properties depend on user choices and no meaningful default value can be set.
// Those properties will be defined later in this file.
// [[codeblock]]
// We enter the namespace Dumux::Properties in order to import the entire Dumux namespace for general use:
namespace Dumux::Properties {

// declaration of the `IncompressibleTest` type tag for the single-phase flow problem
namespace TTag {
struct IncompressibleTest { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
}
// [[/codeblock]]

// ### Property definitions for the type tag `IncompressibleTest`
//
// With the type tag `IncompressibleTest` that we have declared for this single-phase problem,
// we can now define several required compile-time `properties`.
// [[codeblock]]
// We use a structured 2D grid:
template<class TypeTag>
struct Grid<TypeTag, TTag::IncompressibleTest> { using type = Dune::YaspGrid<2>; };

// The problem class specifies initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::IncompressibleTest> { using type = OnePTestProblem<TypeTag>; };

// We define the spatial parameters for our simulation:
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::IncompressibleTest>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

// In the following we define the fluid system to be used (here: liquid water):
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::IncompressibleTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Here, we enable grid-wide caching of geometries and variables. The caches store values
// that were already calculated for later usage. This increases the memory demand but makes the simulation faster.
// This enables grid-wide caching of the volume variables.
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
//This enables grid wide caching for the flux variables.
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
// This enables grid-wide caching for the finite volume grid geometry
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
// [[/codeblock]]

// We use the local residual that contains analytic derivative methods for incompressible flow.
// The one-phase flow model (defined with the `OneP` type tag) uses a default implementation of the
// local residual for single-phase flow. However, in this example we are using an
// incompressible fluid phase. Therefore, we are including the specialized local
// residual which contains functionality to analytically compute the entries of
// the Jacobian matrix for the case of an incompressible fluid phase. We will use this in the main file.
// [[codeblock]]
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::IncompressibleTest> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties
// [[/codeblock]]

#endif
