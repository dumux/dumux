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

#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROPERTIES_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROPERTIES_HH

// ## Compile-time settings (`properties.hh`)
// This file defines the `TypeTag` used for the simulation in this example, for
// which we specialize a number of compile-time `properties`.
// [[content]]
// ### Includes
// [[details]] includes
#include <dune/foamgrid/foamgrid.hh> // for `Dune::FoamGrid`

// The `OneP` type tag specializes most of the `properties` required for single-
// phase flow simulations in DuMu<sup>x</sup>. We will use this in the following to inherit the
// respective properties, and subsequently specialize those properties for our
// type tag, which we want to modify or for which no meaningful default can be set.
#include <dumux/porenetwork/1p/model.hh>// for `TTag::PNMOneP`

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
struct PNMUpscaling { using InheritsFrom = std::tuple<PNMOneP>; };
}

// ### Property specializations
//
// In the following piece of code, mandatory `properties` for which no meaningful
// default can be set, are specialized for our type tag `PNMUpscaling`.
// [[codeblock]]
// We use `dune-foamgrid`, which is especially tailored for 1D networks.
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMUpscaling>
{ using type = Dune::FoamGrid<1, 3>; };

// The problem class specifying initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMUpscaling>
{ using type = UpscalingProblem<TypeTag>; };

//! The spatial parameters to be employed.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMUpscaling>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::UpscalingSpatialParams<GridGeometry, Scalar>;
};

//! The advection type.
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMUpscaling>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, true/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

// We use a single liquid phase consisting of a component with constant fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMUpscaling>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};
// [[/codeblock]]

// Moreover, here we use a local residual specialized for incompressible flow
// that contains functionality related to analytic differentiation.
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMUpscaling>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties
// [[/content]]
#endif
