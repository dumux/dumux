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

// ## Header guard
// The header guard (or include guard) prevents compilation errors due to duplicate definitions.
// Here, a unique name needs to be defined for the header file:
#ifndef DUMUX_LENSPROBLEM_POINTSOURCE_ADAPTIVE_HH
#define DUMUX_LENSPROBLEM_POINTSOURCE_ADAPTIVE_HH

// ## Include files
// The cell centered, two-point-flux discretization scheme is included:
#include <dumux/discretization/cctpfa.hh>

// The fluid properties are specified in the following headers:
#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

// This is the porous medium problem class that this class is derived from:
#include <dumux/porousmediumflow/problem.hh>
// The two-phase flow model is included:
#include <dumux/porousmediumflow/2p/model.hh>
// The local residual for incompressible flow is included:
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

// We include the header that specifies all spatially variable parameters:
#include "spatialparams.hh"

// A container to read values for the initial condition is included:
#include <dumux/io/container.hh>

// ## Define basic properties for our simulation
// We enter the namespace Dumux. All Dumux functions and classes are in a namespace Dumux,
// to make sure they don't clash with symbols from other libraries you may want to use in conjunction with Dumux.
// One could use these functions and classes by prefixing every use of these names by ::,
// but that would quickly become cumbersome and annoying.
// Rather, we simply import the entire Dumux namespace for general use:
namespace Dumux {
  // The problem class is forward declared:
  template<class TypeTag>
  class PointSourceProblem;

  // We enter the namespace Properties, which is a sub-namespace of the namespace Dumux:
  namespace Properties {
  // A TypeTag for our simulation is created which inherits from the two-phase flow model and the
  // cell centered, two-point-flux discretization scheme.
  namespace TTag {
  struct PointSourceExample { using InheritsFrom = std::tuple<TwoP, CCTpfaModel>; };
  }

  // We use non-conforming refinement in our simulation:
  template<class TypeTag>
  struct Grid<TypeTag, TTag::PointSourceExample> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };

  // The problem class specifies initial and boundary conditions:
  template<class TypeTag>
  struct Problem<TypeTag, TTag::PointSourceExample> { using type = PointSourceProblem<TypeTag>; };

  // The local residual contains analytic derivative methods for incompressible flow:
  template<class TypeTag>
  struct LocalResidual<TypeTag, TTag::PointSourceExample> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

  // In the following we define our fluid properties.
  template<class TypeTag>
  struct FluidSystem<TypeTag, TTag::PointSourceExample>
  {
      // We define a convenient shortcut to the property Scalar:
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
      // First, we create a fluid system that consists of one liquid water phase. We use the simple
      // description of water, which means we do not use tabulated values but more general equations of state.
      using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
      // Second, we create another fluid system consisting of a liquid phase as well, the Trichlorethene (DNAPL) phase:
      using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
      // Third, we combine both fluid systems in our final fluid system which consist of two immiscible liquid phases:
      using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
  };

   // we set the formulation for the primary variables to p0s1. In this case that means that the water pressure and the DNAPL saturation are our primary variables.
   template<class TypeTag>
   struct Formulation<TypeTag, TTag::PointSourceExample>
   { static constexpr auto value = TwoPFormulation::p0s1; };

  // We define the spatial parameters for our simulation:
  template<class TypeTag>
  struct SpatialParams<TypeTag, TTag::PointSourceExample>
  {
      // We define convenient shortcuts to the properties FVGridGeometry and Scalar:
  private:
      using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
      // Finally we set the spatial parameters:
  public:
      using type = TwoPTestSpatialParams<FVGridGeometry, Scalar>;
  };

  // We enable caching for the grid volume variables, the grid flux variables and the FV grid geometry. The cache
  // stores values that were already calculated for later usage. This makes the simulation faster.
  template<class TypeTag>
  struct EnableGridVolumeVariablesCache<TypeTag, TTag::PointSourceExample> { static constexpr bool value = false; };
  template<class TypeTag>
  struct EnableGridFluxVariablesCache<TypeTag, TTag::PointSourceExample> { static constexpr bool value = false; };
  template<class TypeTag>
  struct EnableFVGridGeometryCache<TypeTag, TTag::PointSourceExample> { static constexpr bool value = false; };

  //We leave the namespace Properties.
  }

// ## The problem class
// We enter the problem class where all necessary boundary conditions and initial conditions are set for our simulation.
// As this is a porous medium problem, we inherit from the basic PorousMediumFlowProblem.
template <class TypeTag >
class PointSourceProblem : public PorousMediumFlowProblem<TypeTag>
{
    // We use convenient declarations that we derive from the property system.
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimensionworld>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using PointSource =  GetPropType<TypeTag, Properties::PointSource>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    // We define some indices for convenient use in the problem class:
    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationDNAPLIdx = Indices::saturationIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        dnaplPhaseIdx = FluidSystem::phase1Idx
    };

public:
    // This is the constructor of our problem class:
    PointSourceProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
  : ParentType(fvGridGeometry)
  {
    // We read in the values for the initial condition of our simulation:
    initialValues_ = readFileToContainer<std::vector<PrimaryVariables>>("initialsolutioncc.txt");
  }

  // First, we define the type of boundary conditions depending on location. Two types of boundary conditions
  // can be specified: Dirichlet or Neumann boundary condition. On a Dirichlet boundary, the values of the
  // primary variables need to be fixed. On a Neumann boundary condition, values for derivatives need to be fixed.
  // Mixed boundary conditions (different types for different equations on the same boundary) are not accepted.
  BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
  {
    BoundaryTypes values;
    // We specify Dirichlet boundaries on the left and right hand side of our domain:
    if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
        values.setAllDirichlet();
    else
        // The top and bottom of our domain are Neumann boundaries:
        values.setAllNeumann();
    return values;
  }

      // Second, we specify the values for the Dirichlet boundaries, depending on location. As mentioned,
      // we need to fix values of our two primary variables: the water pressure
      // and the Trichlorethene saturation.
      PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
      {
          // To determine the density of water for a given state, we build a fluid state with the given conditions:
          PrimaryVariables values;
          GetPropType<TypeTag, Properties::FluidState> fluidState;
          fluidState.setTemperature(temperature());
          fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
          fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

          // The density is then calculated by the fluid system:
          Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

          // The water phase pressure is the hydrostatic pressure, scaled with a factor:
          Scalar height = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
          Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];
          Scalar alpha = 1 + 1.5/height;
          Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
          Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;

          values[pressureH2OIdx] = 1e5 - factor*densityW*this->spatialParams().gravity(globalPos)[1]*depth;
          // The saturation of the DNAPL Trichlorethene is zero on our Dirichlet boundary:
          values[saturationDNAPLIdx] = 0.0;

          return values;
      }

      // Third, we specify the values for the Neumann boundaries.
      // In our case, we need to specify mass fluxes for our two liquid phases.
      // Inflow is denoted by a negative sign, outflow by a positive sign.
      NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
      {
          // We initialize the fluxes with zero:
          NumEqVector values(0.0);
          // At the inlet, we specify an inflow for our DNAPL Trichlorethene.
          // The units are kg/(m^2 s).
          if (onInlet_(globalPos))
              values[contiDNAPLEqIdx] = -0.04;

          return values;
      }

  // Last, we specify the initial conditions. The initial condition need to be set for all primary variables.
  // Here, we take the data from the file that we read in previously.
  PrimaryVariables initial(const Element& element) const
  {
      // The input data is written for a uniform grid with discretization length delta.
      // Accordingly, we need to find the index of our cells, depending on the x and y coordinates,
      // that corresponds to the indices of the input data set.
      const auto delta = 0.0625;
      unsigned int cellsX = this->fvGridGeometry().bBoxMax()[0]/delta;
      const auto globalPos = element.geometry().center();

      unsigned int dataIdx = std::trunc(globalPos[1]/delta) * cellsX + std::trunc(globalPos[0]/delta);
      return initialValues_[dataIdx];
  }

  // We need to specify a constant temperature for our isothermal problem.
  // Fluid properties that depend on temperature will be calculated with this value.
  Scalar temperature() const
  {
      return 293.15; // 10°C
  }

  // Additionally, we set a point source. The point source can be solution dependent.
  // It is specified in form of a vector that contains source values for alle phases and positions in space.
  // The first entry is a tupel containing the position in space, the second entry contains a tupel with the source (unit kg/s)
  // for the phases (first phase is the water phase, the second phase is the DNAPL Trichlorethene phase).
  void addPointSources(std::vector<PointSource>& pointSources) const
  {
      pointSources.push_back(PointSource({0.502, 3.02}, {0, 0.1}));
  }

  // We define private global functions that are used to determine if a point in space is on the left, right or upper boundary, or
  // at the inlet.
  private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
        Scalar lambda = (this->fvGridGeometry().bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    // Our private global variables are the epsilon value and the vector containing the initial values read from file.
    static constexpr Scalar eps_ = 1e-6;
    std::vector<PrimaryVariables> initialValues_;

    // This is everything the problem class contains.
};
// We leave the namespace Dumux here, too.
}

#endif
