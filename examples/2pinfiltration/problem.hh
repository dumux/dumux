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
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENSPROBLEM_POINTSOURCE_ADAPTIVE_HH
#define DUMUX_LENSPROBLEM_POINTSOURCE_ADAPTIVE_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "spatialparams.hh"

#include <dumux/io/container.hh>

namespace Dumux {
  // forward declarations
  template<class TypeTag> class PointSourceTestProblem;

  namespace Properties {
  // Create new type tags
  namespace TTag {
  struct TwoPAdaptivePointSource { using InheritsFrom = std::tuple<TwoP, CCTpfaModel>; };
  } // end namespace TTag

  //! Use non-conforming refinement
  #if HAVE_DUNE_ALUGRID
  template<class TypeTag>
  struct Grid<TypeTag, TTag::TwoPAdaptivePointSource> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
  #else
  template<class TypeTag>
  struct Grid<TypeTag, TTag::TwoPAdaptivePointSource> { using type = Dune::YaspGrid<2>; };
  #endif

  template<class TypeTag>
  struct Problem<TypeTag, TTag::TwoPAdaptivePointSource> { using type = PointSourceTestProblem<TypeTag>; };

  // the local residual containing the analytic derivative methods
  template<class TypeTag>
  struct LocalResidual<TypeTag, TTag::TwoPAdaptivePointSource> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

  // Set the fluid system
  template<class TypeTag>
  struct FluidSystem<TypeTag, TTag::TwoPAdaptivePointSource>
  {
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
      using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
      using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
      using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
  };

  // Set the spatial parameters
  template<class TypeTag>
  struct SpatialParams<TypeTag, TTag::TwoPAdaptivePointSource>
  {
  private:
      using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
  public:
      using type = TwoPTestSpatialParams<FVGridGeometry, Scalar>;
  };

  // Enable caching
  template<class TypeTag>
  struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPAdaptivePointSource> { static constexpr bool value = false; };
  template<class TypeTag>
  struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPAdaptivePointSource> { static constexpr bool value = false; };
  template<class TypeTag>
  struct EnableFVGridGeometryCache<TypeTag, TTag::TwoPAdaptivePointSource> { static constexpr bool value = false; };

  } // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */
template <class TypeTag >
class PointSourceTestProblem : public PorousMediumFlowProblem<TypeTag>
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
  PointSourceTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
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
          values[pressureH2OIdx] = 1e5 - factor*densityW*this->gravity()[1]*depth;
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

          // in the test with the oil wet lens, use higher injection rate
          if (this->spatialParams().lensIsOilWet())
              values[contiDNAPLEqIdx] *= 10;

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
   * \brief Evaluates the initial values for a control volume.
   *
   * \param globalPos The global position
   */
  PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
  {
      PrimaryVariables values;
      GetPropType<TypeTag, Properties::FluidState> fluidState;
      fluidState.setTemperature(temperature());
      fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
      fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

      Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

      Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];

      // hydrostatic pressure
      values[pressureH2OIdx] = 1e5 - densityW*this->gravity()[1]*depth;
      values[saturationDNAPLIdx] = 0;
      return values;
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

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_;
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
