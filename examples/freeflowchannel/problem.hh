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
// The dune grid interface from YASP grid is included, which is a structured, conforming grid, which can also be used for parallel simulations.
// Next, the properties of the staggered grid (marker-and-cell) discretization scheme are included, which is the spatial discretization used for free-flow simulations in dumux and is summarized in [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/-/blob/master/slides/dumux-course-intro.pdf).
// The single-phase, isothermal Navier-Stokes model and Navier-Stokes problem class that this class is derived from are also included.
// We will simulate the flow of fluid composed of one liquid phase (`1pliquid.hh`), which will have constant fluid properties (density, viscosity,...) (`constant.hh`).
//<details>
//  <summary>Toggle to expand code (file includes):</summary>
//
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>
// </details>
//
// ### Setup basic properties for our simulation
// We setup the DuMux properties for our simulation (click [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/blob/master/slides/dumux-course-properties.pdf) for DuMux course slides on the property system) within the namespace Properties, which is a sub-namespace of Dumux.
// 1. For every test problem, a new `TypeTag` has to be created, which is done within the namespace `TTag` (subnamespace of `Properties`). It inherits from the Navier-Stokes flow model and the staggered-grid discretization scheme.
// 2. The grid is chosen to be a two-dimensional YASP grid.
// 3. We set the `FluidSystem` to be a one-phase liquid with a single component. The class `Component::Constant` refers to a component with constant fluid properties (density, viscosity, ...) that can be set via the input file in the group `[0.Component]` where the number is the identifier given as template argument to the class template `Component::Constant`.
// 4. The problem class `ChannelExampleProblem`, which is forward declared before we enter `namespace Dumux` and defined later in this file, is defined to be the problem used in this test problem (charaterized by the TypeTag `ChannelExample`). The fluid system, which contains information about the properties such as density, viscosity or diffusion coefficient of the fluid we're simulating, is set to a constant one phase liquid.
// 5. We enable caching for the following classes (which stores values that were already calculated for later usage and thus results in higher memory usage but improved CPU speed): the grid volume variables, the grid flux variables, the finite volume grid geometry.
// <details><summary>Toggle to expand code (property definitions):</summary>
//
namespace Dumux {

template <class TypeTag>
class ChannelExampleProblem;

namespace Properties {

namespace TTag {
struct ChannelExample { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
}

template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelExample> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelExample> { using type = Dumux::ChannelExampleProblem<TypeTag> ; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelExample>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
}
// </details>
//
// ### The problem class
// We enter the problem class where all necessary initial and boundary conditions are set for our simulation.
//
// As this is a Stokes problem, we inherit from the basic `NavierStokesProblem`.
// <details><summary>Toggle to expand code:</summary>
//
template <class TypeTag>
class ChannelExampleProblem : public NavierStokesProblem<TypeTag>
{
    // </details>
    //
    // We use convenient declarations that we derive from the property system.
    //<details>
    //  <summary>Toggle to expand code (convenient declarations)</summary>
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
    //  <summary>Toggle to expand code (constructor)</summary>
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
    //  <summary>Toggle to expand code (`boundaryTypesAtPos`)</summary>
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
    //  <summary>Toggle to expand code (`dirichletAtPos`)</summary>
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
    //  <summary>Toggle to expand code (`initialAtPos`)</summary>
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
    //  <summary>Toggle to expand code (`temperature`)</summary>
    //
    Scalar temperature() const
    { return 273.15 + 10; }
private:
    // </details>
    //
    // The inlet is at the left side of the physical domain.
    //<details>
    //  <summary>Toggle to expand code (`isInlet_`)</summary>
    //
    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }
    // </details>
    //
    // The outlet is at the right side of the physical domain.
    //<details>
    //  <summary>Toggle to expand code (`isOutlet_`)</summary>
    //
    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }
    // </details>
    //
    // Finally, private variables are declared:
    //<details>
    //  <summary>Toggle to expand code (private variables)</summary>
    //
    static constexpr Scalar eps_=1e-6;
    Scalar inletVelocity_;
    Scalar outletPressure_;
};
}
#endif
// </details>
//
