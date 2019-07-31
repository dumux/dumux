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

  // In the following we define our fluids.
  template<class TypeTag>
  struct FluidSystem<TypeTag, TTag::PointSourceExample>
  {
      // We define a convenient shortcut to the property Scalar:
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
      // First, we create a fluid system that consists of one liquid phase:
      using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
      // Second, we create another fluid system consisting of a liquid phase:
      using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
      // Third, we combine both fluid systems in our final fluid system which consist of two immiscible phases:
      using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
  };

  // We define the spatial parameters for our simulation.
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

  // Enable caching
  template<class TypeTag>
  struct EnableGridVolumeVariablesCache<TypeTag, TTag::PointSourceExample> { static constexpr bool value = false; };
  template<class TypeTag>
  struct EnableGridFluxVariablesCache<TypeTag, TTag::PointSourceExample> { static constexpr bool value = false; };
  template<class TypeTag>
  struct EnableFVGridGeometryCache<TypeTag, TTag::PointSourceExample> { static constexpr bool value = false; };

  } // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */
template <class TypeTag >
class PointSourceProblem : public PorousMediumFlowProblem<TypeTag>
{
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
    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationDNAPLIdx = Indices::saturationIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        dnaplPhaseIdx = FluidSystem::phase1Idx
    };

public:
    PointSourceProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
  : ParentType(fvGridGeometry)
  {
    initialValues_ = readFileToContainer<std::vector<PrimaryVariables>>("initialsolutioncc.txt");
  }

  /*!
   * \brief Specifies which kind of boundary condition should be
   *        used for which equation on a given boundary segment
   *
   * \param globalPos The global position
   */
  BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
  {
    BoundaryTypes values;
    if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
        values.setAllDirichlet();
    else
        values.setAllNeumann();
    return values;
  }

      /*!
       * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
       *
       * \param globalPos The global position
       */
      PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
      {
          PrimaryVariables values;
          GetPropType<TypeTag, Properties::FluidState> fluidState;
          fluidState.setTemperature(temperature());
          fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
          fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

          Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

          Scalar height = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
          Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];
          Scalar alpha = 1 + 1.5/height;
          Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
          Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;

          // hydrostatic pressure scaled by alpha
          values[pressureH2OIdx] = 1e5 - factor*densityW*this->spatialParams().gravity(globalPos)[1]*depth;
          values[saturationDNAPLIdx] = 0.0;

          return values;
      }

      /*!
       * \brief Evaluates the boundary conditions for a Neumann boundary segment.
       *
       * \param globalPos The position of the integration point of the boundary segment.
       *
       * For this method, the \a values parameter stores the mass flux
       * in normal direction of each phase. Negative values mean influx.
       */
      NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
      {
          NumEqVector values(0.0);
          if (onInlet_(globalPos))
              values[contiDNAPLEqIdx] = -0.04; // kg / (m * s)

          return values;
      }

  /*!
   * \brief Evaluates the initial value for an element for cell-centered models.
   */
  PrimaryVariables initial(const Element& element) const
  {
      const auto delta = 0.0625;
      unsigned int cellsX = this->fvGridGeometry().bBoxMax()[0]/delta;
      const auto globalPos = element.geometry().center();

      // the input data corresponds to a uniform grid with discretization length deltaX_
      unsigned int dataIdx = std::trunc(globalPos[1]/delta) * cellsX + std::trunc(globalPos[0]/delta);
      return initialValues_[dataIdx];
  }

  /*!
   * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
   *
   * This is not specific to the discretization. By default it just
   * throws an exception so it must be overloaded by the problem if
   * no energy equation is used.
   */
  Scalar temperature() const
  {
      return 293.15; // 10°C
  }


    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSources that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in untis
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        // inject 0.1 kg/s of non-wetting phase at position (0.502, 3.02);
        pointSources.push_back(PointSource({0.502, 3.02}, {0, 0.1}));
    }

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

    static constexpr Scalar eps_ = 1e-6;
    std::vector<PrimaryVariables> initialValues_;
};

} // end namespace Dumux

#endif
