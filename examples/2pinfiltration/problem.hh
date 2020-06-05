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

#ifndef DUMUX_EXAMPLE_TWOP_INFILTRATION_PROBLEM_HH
#define DUMUX_EXAMPLE_TWOP_INFILTRATION_PROBLEM_HH

// ## The file `problem.hh`
// [[content]]
//
// ### Includes
// We start with includes for `PorousMediumFlowProblem`, `readFileToContainer`,
// `BoundaryTypes` and `GetPropType` (used below).
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/io/container.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>

// ### Problem class
// The problem class `PointSourceProblem` implements boundary and initial conditions.
// It derives from the `PorousMediumFlowProblem` class.
namespace Dumux {

template <class TypeTag>
class PointSourceProblem : public PorousMediumFlowProblem<TypeTag>
{
    // The class implementation starts with some alias declarations and index definitions for convenience
    // [[details]] local alias declarations and indices
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimensionworld>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using PointSource =  GetPropType<TypeTag, Properties::PointSource>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationDNAPLIdx = Indices::saturationIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        dnaplPhaseIdx = FluidSystem::phase1Idx
    };
    // [[/details]]
    //
    // In the constructor of the class, we call the parent type's constructor
    // and read the intial values for the primary variables from a text file.
    // The function `readFileToContainer` is implemented in the header `dumux/io/container.hh`.
public:
    PointSourceProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        initialValues_ = readFileToContainer<std::vector<PrimaryVariables>>("initialsolutioncc.txt");
    }

    // For isothermal problems, Dumux requires problem classes to implement a `temperature()`
    // member function. Fluid properties that depend on temperature will be calculated with the specified temperature.
    Scalar temperature() const
    {
        return 293.15; // 10°C
    }

    // #### Boundary types
    // We define the type of boundary conditions depending on location. Two types of boundary conditions
    // can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
    // primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.
    // Mixed boundary conditions (different types for different equations on the same boundary) are not accepted.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        // Dirichlet boundaries on the left and right hand side of the domain
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
            values.setAllDirichlet();
        // and Neumann boundaries otherwise (top and bottom of the domain)
        else
            values.setAllNeumann();
        return values;
    }
    // [[/codeblock]]

    // #### Dirichlet boundaries
    // We specify the values for the Dirichlet boundaries, depending on location.
    // We need to fix values for the two primary variables: the water pressure
    // and the DNAPL saturation.
    // [[codeblock]]
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // To determine the density of water for a given state, we build a fluid state with the given conditions:
        PrimaryVariables values;
        GetPropType<TypeTag, Properties::FluidState> fluidState;
        fluidState.setTemperature(temperature());
        fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

        // The density is then calculated by the fluid system:
        const Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

        // The water phase pressure is the hydrostatic pressure, scaled with a factor:
        const Scalar height = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        const Scalar alpha = 1 + 1.5/height;
        const Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        const Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;

        values[pressureH2OIdx] = 1e5 - factor*densityW*this->spatialParams().gravity(globalPos)[1]*depth;
        // The saturation of the DNAPL Trichlorethene is zero on our Dirichlet boundary:
        values[saturationDNAPLIdx] = 0.0;

        return values;
    }
    // [[/codeblock]]

    // #### Neumann boundaries
    // In our case, we need to specify mass fluxes for our two liquid phases.
    // Negative sign means influx and the unit of the boundary flux is $`kg/(m^2 s)`$.
    // On the inlet area, we set a DNAPL influx of $`0.04 kg/(m^2 s)`$. On all other
    // Neumann boundaries, the boundary flux is zero.
    // [[codeblock]]
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos))
            values[contiDNAPLEqIdx] = -0.04;

        return values;
    }
    // [[/codeblock]]

    // #### Initial conditions
    // The initial condition needs to be set for all primary variables.
    // Here, we take the data from the file that we read in previously.
    // [[codeblock]]
    PrimaryVariables initial(const Element& element) const
    {
        // The input data is written for a uniform grid with discretization length delta.
        // Accordingly, we need to find the index of our cells, depending on the x and y coordinates,
        // that corresponds to the indices of the input data set.
        const auto delta = 0.0625;
        const unsigned int cellsX = this->gridGeometry().bBoxMax()[0]/delta;
        const auto globalPos = element.geometry().center();

        const unsigned int dataIdx = std::trunc(globalPos[1]/delta) * cellsX + std::trunc(globalPos[0]/delta);
        return initialValues_[dataIdx];
    }
    // [[/codeblock]]

    // #### Point source
    // In this scenario, we set a point source (e.g. modeling a well). The point source value can be solution dependent.
    // Point sources are added by pushing them into the vector `pointSources`.
    // The `PointSource` constructor takes two arguments.
    // The first argument is a coordinate array containing the position in space,
    // the second argument is an array of source value for each equation (in units of $`kg/s`$).
    // Recall that the first eqution is the water phase mass balance
    // and the second equation is the DNAPL phase mass balance.
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        pointSources.push_back(PointSource({0.502, 3.02}, {0, 0.1}));
    }

    // In the private part of the class, we define some helper functions for
    // the boundary conditions and local variables.
    // [[details]] private data members and functions
private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        Scalar lambda = (this->gridGeometry().bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    static constexpr Scalar eps_ = 1e-6;
    std::vector<PrimaryVariables> initialValues_;
};

} // end namespace Dumux
// [[/details]]
// [[/content]]
#endif
