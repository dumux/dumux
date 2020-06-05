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

#ifndef DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROBLEM_HH
#define DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROBLEM_HH

// ## Initial and boundary conditions (`problem.hh`)
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the Navier-Stokes single-phase flow simulation.
//
// [[content]]
//
// ### Include files
//
// Include the `NavierStokesProblem` class, the base
// class from which we will derive.
#include <dumux/freeflow/navierstokes/problem.hh>
// Include the `NavierStokesBoundaryTypes` class which specifies the boundary types set in this problem.
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As we are solving a problem related to free flow, we inherit from the base class `NavierStokesProblem`.
// [[codeblock]]
namespace Dumux {

template <class TypeTag>
class ChannelExampleProblem : public NavierStokesProblem<TypeTag>
{
    // A few convenience aliases used throughout this class.
    using ParentType = NavierStokesProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<PrimaryVariables::size()>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // This is the constructor of our problem class:
    // Within the constructor, we set the inlet velocity to a run-time specified value.
    // If no run-time value is specified, we set the outlet pressure to 1.1e5 Pa.
    ChannelExampleProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1.1e5);
    }
    // [[/codeblock]]

    // #### Boundary conditions
    // With the following function we define the __type of boundary conditions__ depending on the location.
    // Three types of boundary conditions can be specified: Dirichlet, Neumann or outflow boundary conditions. On
    // Dirichlet boundaries, the values of the primary variables need to be fixed. On a Neumann boundaries,
    // values for derivatives need to be fixed. Outflow conditions set a gradient of zero in normal direction towards the boundary
    // for the respective primary variables (excluding pressure).
    // When Dirichlet conditions are set for the pressure, the velocity gradient
    // with respect to the direction normal to the boundary is automatically set to zero.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if (isInlet_(globalPos))
        {
            // We specify Dirichlet boundary conditions for the velocity on the left of our domain
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else if (isOutlet_(globalPos))
        {
            // We fix the pressure on the right side of the domain
            values.setDirichlet(Indices::pressureIdx);
        }
        else
        {
            // We specify Dirichlet boundary conditions for the velocity on the remaining boundaries (lower and upper wall)
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        return values;
    }
    // [[/codeblock]]

    // The following function specifies the __values on Dirichlet boundaries__.
    // We need to define values for the primary variables (velocity and pressure).
    // [[codeblock]]
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // Use the initial values as default Dirichlet values
        PrimaryVariables values = initialAtPos(globalPos);

        // Set a no-slip condition at the top and bottom wall of the channel
        if (!isInlet_(globalPos))
            values[Indices::velocityXIdx] = 0.0;

        return values;
    }
    // [[/codeblock]]

    // The following function defines the initial conditions.
    // [[codeblock]]
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        // Set the pressure and velocity values
        values[Indices::pressureIdx] = outletPressure_;
        values[Indices::velocityXIdx] = inletVelocity_;
        values[Indices::velocityYIdx] = 0.0;

        return values;
    }
    // [[/codeblock]]

    // #### Temperature distribution
    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    // This would be important if another fluidsystem was used.
    Scalar temperature() const
    { return 273.15 + 10; }

// The inlet is on the left side of the physical domain.
// [[codeblock]]
private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < eps_; }
    // [[/codeblock]]

    // The outlet is on the right side of the physical domain.
    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    // Finally, private variables are declared:
    // [[codeblock]]
    static constexpr Scalar eps_ = 1e-6;
    Scalar inletVelocity_;
    Scalar outletPressure_;

}; // end class definition of ChannelExampleProblem
} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
