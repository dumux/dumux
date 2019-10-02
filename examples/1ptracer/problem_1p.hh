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
#ifndef DUMUX_ONEP_TRACER_TEST_PROBLEM_HH
#define DUMUX_ONEP_TRACER_TEST_PROBLEM_HH

//Before we enter the problem class containing initial and boundary conditions, we include necessary files and introduce properties.
// ### Include files
// The dune grid interphase is included here:
#include <dune/grid/yaspgrid.hh>

// The cell centered, two-point-flux discretization scheme is included:
#include <dumux/discretization/cctpfa.hh>
// The one-phase flow model is included:
#include <dumux/porousmediumflow/1p/model.hh>
// This is the porous medium problem class that this class is derived from:
#include <dumux/porousmediumflow/problem.hh>

// The fluid properties are specified in the following headers:
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// The local residual for incompressible flow is included:
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

// We include the header that specifies all spatially variable parameters:
#include "spatialparams_1p.hh"

// ### Define basic properties for our simulation
// We enter the namespace Dumux in order to import the entire Dumux namespace for general use
namespace Dumux {

// The problem class is forward declared:
template<class TypeTag>
class OnePTestProblem;

// We enter the namespace Properties, which is a sub-namespace of the namespace Dumux:
namespace Properties {
// A `TypeTag` for our simulation is created which inherits from the one-phase flow model and the cell centered, two-point-flux discretization scheme.
namespace TTag {
struct IncompressibleTest { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
}

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
    // We define convenient shortcuts to the properties `FVGridGeometry` and `Scalar`:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    // Finally, we set the spatial parameters:
    using type = OnePTestSpatialParams<FVGridGeometry, Scalar>;
};

// The local residual contains analytic derivative methods for incompressible flow:
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::IncompressibleTest> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

  // In the following we define our fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::IncompressibleTest>
{
    // We define a convenient shortcut to the property Scalar:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    // We create a fluid system that consists of one liquid water phase. We use the simple
    // description of water, which means we do not use tabulated values but more general equations of state.
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// We enable caching for the grid volume variables.
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
// We enable caching for the grid flux variables.
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
// We enable caching for the FV grid geometry
template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::IncompressibleTest> { static constexpr bool value = true; };
//The cache stores values that were already calculated for later usage. This makes the simulation faster.
// We leave the namespace Properties.
}

// ### The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium problem, we inherit from the basic `PorousMediumFlowProblem`.
template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // This is the constructor of our problem class:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    // First, we define the type of boundary conditions depending on location. Two types of boundary  conditions
    // can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
    // primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.
    // Mixed boundary conditions (different types for different equations on the same boundary) are not accepted.
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        // we retreive the global position, i.e. the  vector  including  the  global  coordinates
        // of  the  finite  volume
        const auto globalPos = scvf.ipGlobal();
        // we define a small epsilon value
        Scalar eps = 1.0e-6;
        // We specify Dirichlet boundaries on the top and bottom of our domain:
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->fvGridGeometry().bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
        else
            // The top and bottom of our domain are Neumann boundaries:
            values.setAllNeumann();

        return values;
    }

    // Second, we specify the values for the Dirichlet boundaries. We need to fix values of our  primary variable
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        // we retreive again the global position
        const auto& pos = scvf.ipGlobal();
        PrimaryVariables values(0);
        // we assign pressure values in [Pa] according to a pressure gradient to 1e5 Pa at the top and 1.1e5 Pa at the bottom.
        values[0] = 1.0e+5*(1.1 - pos[dimWorld-1]*0.1);
        return values;
    }

    // We need to specify a constant temperature for our isothermal problem.
    // Fluid properties that depend on temperature will be calculated with this value.
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

    // This is everything the one phase problem class contains.
};

// We leave the namespace Dumux.
} // end namespace Dumux
#endif
