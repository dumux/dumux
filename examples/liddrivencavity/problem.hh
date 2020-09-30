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

#ifndef DUMUX_LIDDRIVENCAVITY_EXAMPLE_PROBLEM_HH
#define DUMUX_LIDDRIVENCAVITY_EXAMPLE_PROBLEM_HH

// ## Initial and boundary conditions (`problem.hh`)
//
// This file contains the __problem class__ which defines the initial and boundary
// conditions for the Navier-Stokes single-phase flow simulation.
//
// [[content]]
//
// ### Include files
//
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

// Include the `NavierStokesProblem` class, the base
// class from which we will derive.
#include <dumux/freeflow/navierstokes/problem.hh>

// Include the `NavierStokesBoundaryTypes` class which specifies the boundary types set in this problem.
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

// ### The problem class
// As we are solving a problem related to free flow, we create a new class called `LidDrivenCavityExampleProblem`
// and let it inherit from the class `NavierStokesProblem`.
// [[codeblock]]
namespace Dumux {
template <class TypeTag>
class LidDrivenCavityExampleProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // Within the constructor, we set the lid velocity to a run-time specified value.
    LidDrivenCavityExampleProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        lidVelocity_ = getParam<Scalar>("Problem.LidVelocity");
    }
    // [[/codeblock]]

    // #### Temperature distribution
    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    // This would be important if another fluidsystem was used.
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

    // #### Boundary conditions
    // With the following function we define the __type of boundary conditions__ depending on the location.
    // Three types of boundary conditions can be specified: Dirichlet, Neumann or outflow boundary conditions. On
    // Dirichlet boundaries, the values of the primary variables need to be fixed. On a Neumann boundaries,
    // values for derivatives need to be fixed. Outflow conditions set a gradient of zero in normal direction towards the boundary
    // for the respective primary variables (excluding pressure).
    // When Dirichlet conditions are set for the pressure, the velocity gradient
    // with respect to the direction normal to the boundary is automatically set to zero.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // We set Dirichlet values for the velocity at each boundary
        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);

        return values;
    }
    // [[/codeblock]]

    // We define a function for setting a fixed Dirichlet pressure value at a given internal cell.
    // This is required for having a defined pressure level in our closed system domain.
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    {
        auto isLowerLeftCell = [&](const SubControlVolume& scv)
        { return scv.dofIndex() == 0; };

        // We set a fixed pressure in one cell
        return (isLowerLeftCell(scv) && pvIdx == Indices::pressureIdx);
    }

    // The following function specifies the __values on Dirichlet boundaries__.
    // We need to define values for the primary variables (velocity and pressure).
    // [[codeblock]]
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = 1.1e+5;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

        // We set the no slip-condition at the top, that means the fluid has the same velocity as the lid
        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            values[Indices::velocityXIdx] = lidVelocity_;

        return values;
    }
    // [[/codeblock]]

    // The following function defines the initial conditions.
    // [[codeblock]]
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = 1.0e+5;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

        return values;
    }
    // [[/codeblock]]
    // the data members of the problem class
    // [[codeblock]]
private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar lidVelocity_;
};

} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
#endif
