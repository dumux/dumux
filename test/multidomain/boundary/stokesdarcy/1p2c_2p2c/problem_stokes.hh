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
 * \ingroup NavierStokesTests
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_STOKES1P2C_SUBPROBLEM_HH
#define DUMUX_STOKES1P2C_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>

namespace Dumux
{
template <class TypeTag>
class StokesSubProblem;

namespace Properties
{
#if !NONISOTHERMAL
NEW_TYPE_TAG(StokesOnePTwoCTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNC));
#else
NEW_TYPE_TAG(StokesOnePTwoCTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNCNI));
#endif


// Set the grid type
SET_TYPE_PROP(StokesOnePTwoCTypeTag, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// The fluid system
SET_PROP(StokesOnePTwoCTypeTag, FluidSystem)
{
  using H2OAir = FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the water phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

SET_INT_PROP(StokesOnePTwoCTypeTag, ReplaceCompEqIdx, 3);

// Use formulation based on mass fractions
SET_BOOL_PROP(StokesOnePTwoCTypeTag, UseMoles, true);

// Set the problem property
SET_TYPE_PROP(StokesOnePTwoCTypeTag, Problem, Dumux::StokesSubProblem<TypeTag> );

SET_BOOL_PROP(StokesOnePTwoCTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(StokesOnePTwoCTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(StokesOnePTwoCTypeTag, EnableGridVolumeVariablesCache, true);

}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

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
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Stokes"), eps_(1e-6), couplingManager_(couplingManager)
    {
        refVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RefVelocity");
        refPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RefPressure");
        refMoleFrac_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.refMoleFrac");
        refTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RefTemperature");

        diffCoeffAvgType_ = StokesDarcyCouplingOptions::stringToEnum(DiffusionCoefficientAveragingType{},
                                                                     getParamFromGroup<std::string>(this->paramGroup(), "Problem.InterfaceDiffusionCoefficientAvg"));
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return refTemperature_; }

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
        BoundaryTypes values;

        const auto& globalPos = scvf.center();

#if NONISOTHERMAL
            values.setNeumann(Indices::energyBalanceIdx);
#endif

        if (onUpperBoundary_(globalPos) || onLeftBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(Indices::conti0EqIdx + 1);
        }

        if (onRightBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(Indices::conti0EqIdx + 1);

#if NONISOTHERMAL
            values.setOutflow(Indices::energyBalanceIdx);
#endif
        }

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::conti0EqIdx + 1);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBJS(Indices::momentumXBalanceIdx);
        }
        return values;
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
        const auto& globalPos = scvf.dofPosition();
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        FluidState fluidState;
        updateFluidStateForBC_(fluidState, elemVolVars[scv].pressure());

        const Scalar density = useMoles ? fluidState.molarDensity(0) : fluidState.density(0);
        const Scalar xVelocity = xVelocity_(globalPos);

        if (onLeftBoundary_(globalPos))
        {
            // rho*v*X at inflow
            values[Indices::conti0EqIdx + 1] = -xVelocity * density * refMoleFrac();
            values[Indices::conti0EqIdx] = -xVelocity * density * (1.0 - refMoleFrac());

#if NONISOTHERMAL
            values[Indices::energyBalanceIdx] = -xVelocity * fluidState.density(0) * fluidState.enthalpy(0);
#endif
        }

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto massFlux = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
            values[Indices::conti0EqIdx] = massFlux[0];
            values[Indices::conti0EqIdx + 1] = massFlux[1];

#if NONISOTHERMAL
            values[Indices::energyBalanceIdx] = couplingManager().couplingData().energyCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
#endif

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
        updateFluidStateForBC_(fluidState, refPressure());

        const Scalar density = FluidSystem::density(fluidState, 0);

        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = refPressure() + density*this->gravity()[1]*(globalPos[1] - this->fvGridGeometry().bBoxMin()[1]);
        values[Indices::conti0EqIdx + 1] = refMoleFrac();
        values[Indices::velocityXIdx] = xVelocity_(globalPos);

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = refTemperature();
#endif

        return values;
    }

    //! \brief Returns the reference velocity.
    const Scalar refVelocity() const
    { return refVelocity_ ;}

    //! \brief Returns the reference pressure.
    const Scalar refPressure() const
    { return refPressure_; }

    //! \brief Returns the reference mass fraction.
    const Scalar refMoleFrac() const
    { return refMoleFrac_; }

    //! \brief Returns the reference temperature.
    const Scalar refTemperature() const
    { return refTemperature_; }


    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar time() const
    { return timeLoop_->time(); }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(scvf);
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
    void updateFluidStateForBC_(FluidState& fluidState, const Scalar pressure) const
    {
        fluidState.setTemperature(refTemperature());
        fluidState.setPressure(0, pressure);
        fluidState.setSaturation(0, 1.0);
        fluidState.setMoleFraction(0, 1, refMoleFrac());
        fluidState.setMoleFraction(0, 0, 1.0 - refMoleFrac());

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, 0);

        const Scalar density = FluidSystem::density(fluidState, paramCache, 0);
        fluidState.setDensity(0, density);

        const Scalar molarDensity = FluidSystem::molarDensity(fluidState, paramCache, 0);
        fluidState.setMolarDensity(0, molarDensity);

        const Scalar enthalpy = FluidSystem::enthalpy(fluidState, paramCache, 0);
        fluidState.setEnthalpy(0, enthalpy);
    }

    //! \brief set the profile of the inflow velocity (horizontal direction)
    const Scalar xVelocity_(const GlobalPosition &globalPos) const
    {
        const Scalar vmax = refVelocity();
        return  4 * vmax * (globalPos[1] - this->fvGridGeometry().bBoxMin()[1]) * (this->fvGridGeometry().bBoxMax()[1] - globalPos[1])
                / (height_() * height_());
    }

    // the height of the free-flow domain
    const Scalar height_() const
    { return this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]; }

    Scalar eps_;

    Scalar refVelocity_;
    Scalar refPressure_;
    Scalar refMoleFrac_;
    Scalar refTemperature_;
    std::string problemName_;
    TimeLoopPtr timeLoop_;

    std::shared_ptr<CouplingManager> couplingManager_;

    DiffusionCoefficientAveragingType diffCoeffAvgType_;
};
} //end namespace

#endif // DUMUX_STOKES1P2C_SUBPROBLEM_HH
