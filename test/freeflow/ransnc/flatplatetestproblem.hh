// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup RANSNCTests
 * \brief Flat plate test for the multi-component staggered grid Reynolds-averaged Navier-Stokes model
 */
#ifndef DUMUX_RANS_NC_TEST_PROBLEM_HH
#define DUMUX_RANS_NC_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

#if LOWREKEPSILON
#include <dumux/freeflow/compositional/lowrekepsilonncmodel.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#elif KEPSILON
#include <dumux/freeflow/compositional/kepsilonncmodel.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>
#elif KOMEGA
#include <dumux/freeflow/compositional/komegancmodel.hh>
#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#elif ONEEQ
#include <dumux/freeflow/compositional/oneeqncmodel.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#else
#include <dumux/freeflow/compositional/zeroeqncmodel.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
#endif

namespace Dumux
{
template <class TypeTag>
class FlatPlateNCTestProblem;

namespace Properties
{

#if NONISOTHERMAL
  #if LOWREKEPSILON
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, LowReKEpsilonNCNI));
  #elif KEPSILON
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, KEpsilonNCNI));
  #elif KOMEGA
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, KOmegaNCNI));
  #elif ONEEQ
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNCNI));
  #else
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNCNI));
  #endif
#else
  #if LOWREKEPSILON
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, LowReKEpsilonNC));
  #elif KEPSILON
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, KEpsilonNC));
  #elif KOMEGA
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, KOmegaNC));
  #elif ONEEQ
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNC));
  #else
  NEW_TYPE_TAG(FlatPlateNCTest, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNC));
  #endif
#endif

// The fluid system
SET_PROP(FlatPlateNCTest, FluidSystem)
{
  using H2OAir = FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the air phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// replace the main component balance eq with a total balance eq
SET_INT_PROP(FlatPlateNCTest, ReplaceCompEqIdx, 0);

// Set the grid type
SET_TYPE_PROP(FlatPlateNCTest, Grid,
              Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(FlatPlateNCTest, Problem, Dumux::FlatPlateNCTestProblem<TypeTag> );

SET_BOOL_PROP(FlatPlateNCTest, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(FlatPlateNCTest, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(FlatPlateNCTest, EnableGridVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(FlatPlateNCTest, UseMoles, true);
} // end namespace Properties

/*!
 * \ingroup RANSNCTests
 * \brief  Test problem for the one-phase model.
 *
 * Dry air is entering from the left side and flows above a 1-D a flat plate.
 * In the middle of the inlet, water vapor is injected, which spreads by turbulent diffusion.
 * For the nonisothermal model the bottom has a constant temperature
 * which is \f$ \unit[30]{K} \f$ higher than the initial and inlet temperature.
 */
template <class TypeTag>
#if LOWREKEPSILON
class FlatPlateNCTestProblem : public LowReKEpsilonProblem<TypeTag>
{
    using ParentType = LowReKEpsilonProblem<TypeTag>;
#elif KEPSILON
class FlatPlateNCTestProblem : public KEpsilonProblem<TypeTag>
{
    using ParentType = KEpsilonProblem<TypeTag>;
#elif KOMEGA
class FlatPlateNCTestProblem : public KOmegaProblem<TypeTag>
{
    using ParentType = KOmegaProblem<TypeTag>;
#elif ONEEQ
class FlatPlateNCTestProblem : public OneEqProblem<TypeTag>
{
    using ParentType = OneEqProblem<TypeTag>;
#else
class FlatPlateNCTestProblem : public ZeroEqProblem<TypeTag>
{
    using ParentType = ZeroEqProblem<TypeTag>;
#endif

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto dimWorld = GET_PROP_TYPE(TypeTag, GridView)::dimensionworld;
    static constexpr auto transportEqIdx = Indices::conti0EqIdx + 1;
    static constexpr auto transportCompIdx = Indices::conti0EqIdx + 1;

public:
    FlatPlateNCTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");

        FluidSystem::init();
        Dumux::TurbulenceProperties<Scalar, dimWorld, true> turbulenceProperties;
        FluidState fluidState;
        const auto phaseIdx = 0;
        fluidState.setPressure(phaseIdx, 1e5);
        fluidState.setTemperature(temperature());
        fluidState.setMassFraction(phaseIdx, phaseIdx, 1.0);
        Scalar density = FluidSystem::density(fluidState, phaseIdx);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, phaseIdx) / density;
        Scalar diameter = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        viscosityTilde_ = 1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, diameter, kinematicViscosity);
        turbulentKineticEnergy_ = turbulenceProperties.turbulentKineticEnergy(inletVelocity_, diameter, kinematicViscosity);
#if KOMEGA
        dissipation_ = turbulenceProperties.dissipationRate(inletVelocity_, diameter, kinematicViscosity);
#else
        dissipation_ = turbulenceProperties.dissipation(inletVelocity_, diameter, kinematicViscosity);
#endif
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteRestartFile() const
    {
        return false;
    }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * The isothermal problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(isInlet_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setDirichlet(transportCompIdx);

#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif

#if KEPSILON || KOMEGA || LOWREKEPSILON
            values.setDirichlet(Indices::turbulentKineticEnergyIdx);
            values.setDirichlet(Indices::dissipationIdx);
#endif
#if ONEEQ
            values.setDirichlet(Indices::viscosityTildeIdx);
#endif
        }
        else if(isOutlet_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(transportEqIdx);

#if NONISOTHERMAL
            values.setOutflow(Indices::energyBalanceIdx);
#endif

#if KEPSILON || KOMEGA || LOWREKEPSILON
            values.setOutflow(Indices::turbulentKineticEnergyEqIdx);
            values.setOutflow(Indices::dissipationEqIdx);
#endif
#if ONEEQ
            values.setOutflow(Indices::viscosityTildeIdx);
#endif
        }
        else if(isOnWallAtPos(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(transportEqIdx);

#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif

#if KEPSILON || LOWREKEPSILON
            values.setDirichlet(Indices::turbulentKineticEnergyEqIdx);
            values.setDirichlet(Indices::dissipationEqIdx);
#elif KOMEGA
            values.setDirichlet(Indices::turbulentKineticEnergyEqIdx);
#endif
#if ONEEQ
            values.setDirichlet(Indices::viscosityTildeIdx);
#endif
        }
        else
        {
            values.setAllSymmetry();
        }

        return values;
    }

#if KOMEGA
    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     */
    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell(const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const SubControlVolume& scv,
                         int pvIdx) const
    {
        // set a fixed dissipation (omega) for all cells at the wall
        for (const auto& scvf : scvfs(fvGeometry))
            if (isOnWallAtPos(scvf.center()) && pvIdx == Indices::dissipationIdx)
                return true;

        return false;
    }
#endif

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet values at the boundary.
     *
     * \param element The finite element
     * \param scvf the sub control volume face
     * \note used for cell-centered discretization schemes
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        const auto globalPos = scvf.ipGlobal();
        PrimaryVariables values = initialAtPos(globalPos);

        if (time() > 10.0)
        {
            if (isInlet_(globalPos)
                && globalPos[1] > 0.4 * this->fvGridGeometry().bBoxMax()[1]
                && globalPos[1] < 0.6 * this->fvGridGeometry().bBoxMax()[1])
            {
                values[transportCompIdx] = 1e-3;
            }
#if NONISOTHERMAL
            if (isOnWallAtPos(globalPos))
            {
                values[Indices::temperatureIdx] += 30.;
            }
#endif
        }

        return values;
    }

     /*!
      * \brief Evaluate the boundary conditions for fixed values at cell centers
      *
      * \param element The finite element
      * \param scv the sub control volume
      * \note used for cell-centered discretization schemes
      */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolume &scv) const
    {
        const auto globalPos = scv.center();
        PrimaryVariables values(initialAtPos(globalPos));
#if KOMEGA
        using std::pow;
        unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto wallDistance = ParentType::wallDistance_[elementIdx];
        values[Indices::dissipationEqIdx] = 6.0 * ParentType::kinematicViscosity_[elementIdx]
                                            / (ParentType::betaOmega() * pow(wallDistance, 2));
#endif
        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1.0e+5;
        values[transportCompIdx] = 0.0;
#if NONISOTHERMAL
        values[Indices::temperatureIdx] = temperature();
#endif

        // block velocity profile
        values[Indices::velocityXIdx] = 0.0;
        if (!isOnWallAtPos(globalPos))
            values[Indices::velocityXIdx] =  inletVelocity_;
        values[Indices::velocityYIdx] = 0.0;

#if KEPSILON || KOMEGA || LOWREKEPSILON
        values[Indices::turbulentKineticEnergyEqIdx] = turbulentKineticEnergy_;
        values[Indices::dissipationEqIdx] = dissipation_;
        if (isOnWallAtPos(globalPos))
        {
            values[Indices::turbulentKineticEnergyEqIdx] = 0.0;
            values[Indices::dissipationEqIdx] = 0.0;
        }
#endif

#if ONEEQ
        values[Indices::viscosityTildeIdx] = viscosityTilde_;
        if (isOnWallAtPos(globalPos))
        {
            values[Indices::viscosityTildeIdx] = 0.0;
        }
#endif

        return values;
    }

    // \}

    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

    bool isOnWallAtPos(const GlobalPosition& globalPos) const
    {
        return globalPos[1] < eps_;
    }

private:


    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    const Scalar eps_;
    Scalar inletVelocity_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    TimeLoopPtr timeLoop_;
};
} //end namespace

#endif
