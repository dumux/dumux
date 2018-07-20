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
  * \brief The free-flow sub problem
  */
#ifndef DUMUX_RANS1P2C_SUBPROBLEM_HH
#define DUMUX_RANS1P2C_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

#if ONEEQ
#include <dumux/freeflow/compositional/oneeqncmodel.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#elif KEPSILON
#include <dumux/freeflow/compositional/oneeqncmodel.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#elif KOMEGA
#include <dumux/freeflow/compositional/oneeqncmodel.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#elif LOWREKEPSILON
#include <dumux/freeflow/compositional/oneeqncmodel.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#else
#include <dumux/freeflow/compositional/zeroeqncmodel.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
#endif

namespace Dumux
{
template <class TypeTag>
class FreeFlowSubProblem;

namespace Properties
{
#if ONEEQ
NEW_TYPE_TAG(RANSTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNCNI));
#elif KEPSILON
NEW_TYPE_TAG(RANSTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNCNI));
#elif KOMEGA
NEW_TYPE_TAG(RANSTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNCNI));
#elif LOWREKEPSILON
NEW_TYPE_TAG(RANSTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, OneEqNCNI));
#else
NEW_TYPE_TAG(RANSTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEqNCNI));
#endif

// Set the grid type
SET_TYPE_PROP(RANSTypeTag, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// The fluid system
SET_PROP(RANSTypeTag, FluidSystem)
{
  using H2OAir = FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the air phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

SET_INT_PROP(RANSTypeTag, ReplaceCompEqIdx, 3);

// Use formulation based on mass fractions
SET_BOOL_PROP(RANSTypeTag, UseMoles, true);

// Set the problem property
SET_TYPE_PROP(RANSTypeTag, Problem, Dumux::FreeFlowSubProblem<TypeTag> );

SET_BOOL_PROP(RANSTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(RANSTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(RANSTypeTag, EnableGridVolumeVariablesCache, true);

SET_BOOL_PROP(RANSTypeTag, EnableInertiaTerms, true);
}

/*!
 * \brief The free-flow sub problem
 */
template <class TypeTag>
#if ONEEQ
class FreeFlowSubProblem : public OneEqProblem<TypeTag>
{
    using ParentType = OneEqProblem<TypeTag>;
#elif KEPSILON
class FreeFlowSubProblem : public OneEqProblem<TypeTag>
{
    using ParentType = OneEqProblem<TypeTag>;
#elif KOMEGA
class FreeFlowSubProblem : public OneEqProblem<TypeTag>
{
    using ParentType = OneEqProblem<TypeTag>;
#elif LOWREKEPSILON
class FreeFlowSubProblem : public OneEqProblem<TypeTag>
{
    using ParentType = OneEqProblem<TypeTag>;
#else
class FreeFlowSubProblem : public ZeroEqProblem<TypeTag>
{
    using ParentType = ZeroEqProblem<TypeTag>;
#endif

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    static constexpr bool useMoles = GET_PROP_TYPE(TypeTag, ModelTraits)::useMoles();

public:
    FreeFlowSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "RANS"), eps_(1e-6), couplingManager_(couplingManager)
    {
        inletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletVelocity");
        inletPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletPressure");
        inletMoleFrac_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletMoleFrac");
        inletTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletTemperature");

        diffCoeffAvgType_ = StokesDarcyCouplingOptions::stringToEnum(DiffusionCoefficientAveragingType{},
                                                                     getParamFromGroup<std::string>(this->paramGroup(), "Problem.InterfaceDiffusionCoefficientAvg"));

        Dumux::TurbulenceProperties<Scalar, FVGridGeometry::GridView::dimensionworld, true> turbulenceProperties;
        FluidState fluidState;
        fluidState.setPressure(0, inletPressure_);
        fluidState.setMoleFraction(0, 1, inletMoleFrac_);
        fluidState.setMoleFraction(0, 0, 1.0 - inletMoleFrac_);
        fluidState.setTemperature(inletTemperature_);
        Scalar density = FluidSystem::density(fluidState, 0);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0) / density;
        Scalar radius = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        // ideally the viscosityTilde parameter as inflow for the Spalart-Allmaras model should be zero
        viscosityTilde_ = getParam<Scalar>("Problem.InletViscosityTilde",
                                           1e-3 * turbulenceProperties.viscosityTilde(inletVelocity_, radius*2.0, kinematicViscosity, true));
        turbulentKineticEnergy_ = getParam<Scalar>("Problem.InletTurbulentKineticEnergy",
                                                   turbulenceProperties.turbulentKineticEnergy(inletVelocity_, radius*2.0, kinematicViscosity, true));
#if KOMEGA
        dissipation_ = getParam<Scalar>("Problem.InletDissipationRate",
                                        turbulenceProperties.dissipationRate(inletVelocity_, radius*2.0, kinematicViscosity, true));
#else
        dissipation_ = getParam<Scalar>("Problem.InletDissipation",
                                        turbulenceProperties.dissipation(inletVelocity_, radius*2.0, kinematicViscosity, true));
#endif
        std::cout << std::endl;
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return inletTemperature_; }

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes bTypes;

        const auto& globalPos = scvf.center();

#if ONEEQ
            bTypes.setDirichlet(Indices::viscosityTildeIdx);
#endif

        if (onLeftBoundary_(globalPos))
        {
            bTypes.setDirichlet(Indices::velocityXIdx);
            bTypes.setDirichlet(Indices::velocityYIdx);
            bTypes.setDirichlet(Indices::conti0EqIdx + 1);
            bTypes.setDirichlet(Indices::energyBalanceIdx);
        }

        if (onLowerBoundary_(globalPos))
        {
            bTypes.setDirichlet(Indices::velocityXIdx);
            bTypes.setDirichlet(Indices::velocityYIdx);
            bTypes.setNeumann(Indices::conti0EqIdx);
            bTypes.setNeumann(Indices::conti0EqIdx + 1);
            bTypes.setNeumann(Indices::energyBalanceIdx);
        }

        if (onUpperBoundary_(globalPos))
        {
            bTypes.setAllSymmetry();
        }

        if (onRightBoundary_(globalPos))
        {
            bTypes.setDirichlet(Indices::pressureIdx);
            bTypes.setOutflow(Indices::conti0EqIdx + 1);
            bTypes.setOutflow(Indices::energyBalanceIdx);
#if ONEEQ
            bTypes.setOutflow(Indices::viscosityTildeIdx);
#endif
        }

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            bTypes.setCouplingNeumann(Indices::conti0EqIdx);
            bTypes.setCouplingNeumann(Indices::conti0EqIdx + 1);
            bTypes.setCouplingNeumann(Indices::energyBalanceIdx);
            bTypes.setCouplingNeumann(Indices::momentumYBalanceIdx);
            bTypes.setBJS(Indices::momentumXBalanceIdx);
        }
        return bTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element
     * \param scvf The subcontrolvolume face
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(pos);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto massFlux = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
            values[Indices::conti0EqIdx] = massFlux[0];
            values[Indices::conti0EqIdx + 1] = massFlux[1];
            values[Indices::energyBalanceIdx] = couplingManager().couplingData().energyCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
        }
        return values;
    }

    // \}

    /*!
     * \brief Set the coupling manager
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    /*!
     * \brief Get the coupling manager
     */
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    bool isOnWall(const GlobalPosition& globalPos) const
    {
        return (onLowerBoundary_(globalPos));
    }

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
        FluidState fluidState;
        updateFluidStateForBC_(fluidState, inletTemperature(), inletPressure(), inletMoleFrac());

        const Scalar density = FluidSystem::density(fluidState, 0);

        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = inletPressure() + density*this->gravity()[1]*(globalPos[1] - this->fvGridGeometry().bBoxMin()[1]);
        values[Indices::conti0EqIdx + 1] = inletMoleFrac();
        values[Indices::velocityXIdx] = inletVelocity();
        values[Indices::temperatureIdx] = inletTemperature();

        if(onLowerBoundary_(globalPos))
            values[Indices::velocityXIdx] = 0.0;

#if ONEEQ
        values[Indices::viscosityTildeIdx] = viscosityTilde_;
        if (isOnWall(globalPos))
        {
            values[Indices::viscosityTildeIdx] = 0.0;
        }
#endif

        return values;
    }

    //! \brief Returns the inlet velocity.
    const Scalar inletVelocity() const
    { return inletVelocity_ ;}

    //! \brief Returns the inlet pressure.
    const Scalar inletPressure() const
    { return inletPressure_; }

    //! \brief Returns the inlet mass fraction.
    const Scalar inletMoleFrac() const
    { return inletMoleFrac_; }

    //! \brief Returns the inlet temperature.
    const Scalar inletTemperature() const
    { return inletTemperature_; }


    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().permeabilityAtPos(scvf.center());
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
    }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    //! \brief updates the fluid state to obtain required quantities for IC/BC
    void updateFluidStateForBC_(FluidState& fluidState, const Scalar temperature,
                                const Scalar pressure, const Scalar moleFraction) const
    {
        fluidState.setTemperature(temperature);
        fluidState.setPressure(0, pressure);
        fluidState.setSaturation(0, 1.0);
        fluidState.setMoleFraction(0, 1, moleFraction);
        fluidState.setMoleFraction(0, 0, 1.0 - moleFraction);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, 0);

        const Scalar density = FluidSystem::density(fluidState, paramCache, 0);
        fluidState.setDensity(0, density);

        const Scalar molarDensity = FluidSystem::molarDensity(fluidState, paramCache, 0);
        fluidState.setMolarDensity(0, molarDensity);

        const Scalar enthalpy = FluidSystem::enthalpy(fluidState, paramCache, 0);
        fluidState.setEnthalpy(0, enthalpy);
    }

    // the height of the free-flow domain
    const Scalar height_() const
    { return this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]; }

    Scalar eps_;

    Scalar inletVelocity_;
    Scalar inletPressure_;
    Scalar inletMoleFrac_;
    Scalar inletTemperature_;
    Scalar viscosityTilde_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;

    TimeLoopPtr timeLoop_;

    std::shared_ptr<CouplingManager> couplingManager_;

    DiffusionCoefficientAveragingType diffCoeffAvgType_;
};
} //end namespace

#endif // DUMUX_STOKES1P2C_SUBPROBLEM_HH
