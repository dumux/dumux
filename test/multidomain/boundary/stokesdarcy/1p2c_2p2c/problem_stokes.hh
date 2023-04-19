// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_STOKES1P2C_SUBPROBLEM_HH
#define DUMUX_STOKES1P2C_SUBPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingdata.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFaceVariables = typename GetPropType<TypeTag, Properties::GridFaceVariables>::LocalView;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    static constexpr bool useMoles = GetPropType<TypeTag, Properties::ModelTraits>::useMoles();

public:
    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Stokes"), eps_(1e-6), couplingManager_(couplingManager)
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
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    /*!
     * \brief Returns the sources within the domain.
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
            values.setNeumann(Indices::energyEqIdx);
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
            values.setOutflow(Indices::energyEqIdx);
#endif
        }

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::conti0EqIdx + 1);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBeaversJoseph(Indices::momentumXBalanceIdx);
        }
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(pos);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
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
            values[Indices::energyEqIdx] = -xVelocity * fluidState.density(0) * fluidState.enthalpy(0);
#endif
        }

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto massFlux = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
            values[Indices::conti0EqIdx] = massFlux[0];
            values[Indices::conti0EqIdx + 1] = massFlux[1];

#if NONISOTHERMAL
            values[Indices::energyEqIdx] = couplingManager().couplingData().energyCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
#endif

        }
        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        FluidState fluidState;
        updateFluidStateForBC_(fluidState, refPressure());

        const Scalar density = FluidSystem::density(fluidState, 0);

        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = refPressure() + density*this->gravity()[1]*(globalPos[1] - this->gridGeometry().bBoxMin()[1]);
        values[Indices::conti0EqIdx + 1] = refMoleFrac();
        values[Indices::velocityXIdx] = xVelocity_(globalPos);

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = refTemperature();
#endif

        return values;
    }

    //! Returns the reference velocity.
    const Scalar refVelocity() const
    { return refVelocity_ ;}

    //! Returns the reference pressure.
    const Scalar refPressure() const
    { return refPressure_; }

    //! Returns the reference mass fraction.
    const Scalar refMoleFrac() const
    { return refMoleFrac_; }

    //! Returns the reference temperature.
    const Scalar refTemperature() const
    { return refTemperature_; }


    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar time() const
    { return timeLoop_->time(); }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
     *        for the Beavers-Joseph-Saffman boundary condition.
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element, scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the
     *        Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
    }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    //! Updates the fluid state to obtain required quantities for IC/BC
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

    //! Set the profile of the inflow velocity (horizontal direction).
    const Scalar xVelocity_(const GlobalPosition &globalPos) const
    {
        const Scalar vmax = refVelocity();
        return  4 * vmax * (globalPos[1] - this->gridGeometry().bBoxMin()[1]) * (this->gridGeometry().bBoxMax()[1] - globalPos[1])
                / (height_() * height_());
    }

    // the height of the free-flow domain
    const Scalar height_() const
    { return this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]; }

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
} // end namespace Dumux

#endif
