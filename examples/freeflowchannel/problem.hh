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

// ### Header guard
#ifndef DUMUX_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_TEST_PROBLEM_HH

//Before we enter the problem class containing initial and boundary conditions, we include necessary files and introduce properties.
// ### Include files
// The dune grid interface (from YASP) is included, as are the staggered grid discretization scheme, the freeflow model
// and the freeflow Navier-Stokes problem class that this class is derived from. The material (fluid) properties are specified.
//<details>
//  <summary>Click to toggle details</summary>
//
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
// </details>
//
// ### Define basic properties for our simulation
// Basis properties of the simulation are defined, e.g. the model, discretization scheme, grid, fluid properties, caching.
//<details>
//  <summary>Click to toggle details</summary>
//
// We enter the namespace Dumux in order to import the entire Dumux namespace for general use
namespace Dumux {

// The problem class is forward declared:
template <class TypeTag>
class ChannelExampleProblem;

// We enter the namespace Properties, which is a sub-namespace of the namespace Dumux:
namespace Properties {
// Create new type tags
namespace TTag {
// A `TypeTag` for our simulation is created which inherits from the Navier-Stokes flow model and the staggered-grid discretization scheme.
struct ChannelExample { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// We use a structured 2D grid:
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelExample> { using type = Dune::YaspGrid<2>; };

// The problem class specifies initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelExample> { using type = Dumux::ChannelExampleProblem<TypeTag> ; };

// This is where we define the fluid system, which contains information about the properties of the fluid we're simulating. To define the fluid system we first define the property Scalar. We then use this type to create a fluid system that consists of an incompressible fluid of constant visosity.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelExample>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// We enable caching for the grid volume variables.
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
// We enable caching for the grid flux variables.
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
// We enable caching for the FV grid geometry
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
//The cache stores values that were already calculated for later usage. This makes the simulation faster.
// We leave the namespace Properties.
}
// </details>
//
// ### The problem class
// We enter the problem class where all necessary initial and boundary conditions are set for our simulation.
//
//<details>
//  <summary>Click to toggle details</summary>
//
// As this is a Stokes problem, we inherit from the basic `NavierStokesProblem`.
template <class TypeTag>
class ChannelExampleProblem : public NavierStokesProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    using ParentType = NavierStokesProblem<TypeTag>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // </details>
    //
    // There follows the constructor of our problem class:
    // Within the constructor, we set the inlet velocity to a run-time specified value.
    // As no run-time value is specified, we set the outlet pressure to 1.1e5 Pa.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    ChannelExampleProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.1e5);
    }
    // </details>
    //
    // Now, we define the type of initial and boundary conditions depending on location.
    // Two types of boundary  conditions can be specified: Dirichlet and Neumann. On a Dirichlet boundary,
    // the values of the primary variables need to be fixed.
    // On a Neumann boundary condition, values for derivatives need to be fixed.
    // When Dirichlet conditions are set for the pressure, the derivative of the velocity
    // vector with respect to the  direction normal to the boundary is automatically set to
    // zero. This boundary condition is called in-/outflow boundary condition in Dumux.
    // In the following we specify Dirichlet boundaries for velocity on the left of our domain
    // if isInlet_ is true, Dirichlet boundaries for pressure on the right of our domain
    // if isOutlet_ is true and specify Dirichlet boundaries for velocity on the top and bottom
    // of our domain else.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(isInlet_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else if(isOutlet_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
        }
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        return values;
    }
    // </details>
    //
    // Second, we specify the values for the Dirichlet boundaries. We need to fix the values of our primary variables.
    // To ensure a no-slip boundary condition at the top and bottom of the channel, the Dirichlet velocity
    // in x-direction is set to zero if not at the inlet.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        if(!isInlet_(globalPos))
        {
            values[Indices::velocityXIdx] = 0.0;
        }

        return values;
    }
    // </details>
    //
    // We specify the values for the initial conditions.
    // We assign constant values for pressure and velocity components.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        values[Indices::pressureIdx] = outletPressure_;
        values[Indices::velocityXIdx] = inletVelocity_;
        values[Indices::velocityYIdx] = 0.0;

        return values;
    }
    // </details>
    //
    // We need to specify a constant temperature for our isothermal problem.
    // We set it to 10°C.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    Scalar temperature() const
    { return 273.15 + 10; }
private:
    // </details>
    //
    // The inlet is at the left side of the physical domain.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }
    // </details>
    //
    // The outlet is at the right side of the physical domain.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }
    // </details>
    //
    // Finally, private variables are declared:
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    static constexpr Scalar eps_=1e-6;
    Scalar inletVelocity_;
    Scalar outletPressure_;
    // </details>
    //
    // This is everything the freeflow channel problem class contains.
};
// We leave the namespace Dumux.
} // end namespace Dumux
#endif
// </details>
//
