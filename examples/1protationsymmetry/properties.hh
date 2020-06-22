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

#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_PROPERTIES_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_PROPERTIES_HH

// ## Compile-time settings (`properties.hh`)
// This file defines the `TypeTag` used for the simulation in this example, for
// which we specialize a number of compile-time `properties`.
// [[content]]
// ### Includes
// [[details]] includes
#include <dune/grid/yaspgrid.hh> // for `Dune::YaspGrid`
#include <dumux/discretization/box.hh> // for `TTag::BoxModel`


// The `OneP` type tag specializes most of the `properties` required for single-
// phase flow simulations in DuMu<sup>x</sup>. We will use this in the following to inherit the
// respective properties, and subsequently specialize those properties for our
// type tag, which we want to modify or for which no meaningful default can be set.
#include <dumux/porousmediumflow/1p/model.hh> // for `TTag::OneP`

// The local residual for incompressible flow is included.
// The one-phase flow model (included above) uses a default implementation of the
// local residual for single-phase flow. However, in this example we are using an
// incompressible fluid phase. Therefore, we are including the specialized local
// residual which contains functionality to analytically compute the entries of
// the Jacobian matrix. We will use this in the main file.
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

// We will use a single liquid phase consisting of a component with constant fluid properties.
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// As mentioned at the beginning of the documentation of this example, DuMu<sup>x</sup>
// provides specialized implementations of the extrusion type for rotational symmetric geometries.
// The extrusion class takes care of adjusting volume and area
// computations to account for the extrusion about the symmetry axes.
#include <dumux/discretization/extrusion.hh>

// The classes that define the problem and parameters used in this simulation
#include "problem.hh"
#include "spatialparams.hh"
// [[/details]]
//
// ### `TypeTag` definition
// A `TypeTag` for our simulation is defined, which inherits properties from the
// single-phase flow model and the box scheme.
namespace Dumux::Properties {
namespace TTag {
struct OnePRotSym { using InheritsFrom = std::tuple<OneP, BoxModel>; };
}

// ### Property specializations
//
// In the following piece of code, mandatory `properties` for which no meaningful
// default can be set, are specialized for our type tag `OnePRotSym`.
// [[codeblock]]
// We use a structured 1D grid with an offset. This allows us to define the
// computational domain to be between the radii $`r_1`$ and $`r_2`$ as illustrated
// in the beginning of the documentation of this example
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePRotSym>
{ using type =  Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>; };

// The problem class specifying initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePRotSym>
{ using type = RotSymExampleProblem<TypeTag>; };

// Our spatial parameters class defining the permeability and porosity of the porous medium:
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePRotSym>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = RotSymExampleSpatialParams<GridGeometry, Scalar>;
};

// We use a single liquid phase consisting of a component with constant fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePRotSym>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};
// [[/codeblock]]

// As mentioned before, we modify the areas and volumes used
// for integrals over the control volume and the control volume faces by changing the extrusion type.
// Here, we pass these traits to the grid geometry of the box scheme (the scheme
// that we use here) and specialize the `GridGeometry` property accordingly.
// [[codeblock]]
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePRotSym>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    // We take the default traits as basis and exchange the extrusion type
    // The first axis (x-axis) is the radial axis, hence the zero. That means we rotate about the second axis (y-axis).
    struct GGTraits : public BoxDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

public:
    // Pass the above traits to the box grid geometry such that it uses the
    // rotation-symmetric sub-control volumes and faces.
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};
// [[/codeblock]]

// Moreover, here we use a local residual specialized for incompressible flow
// that contains functionality related to analytic differentiation.
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePRotSym>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties
// [[/content]]
#endif
