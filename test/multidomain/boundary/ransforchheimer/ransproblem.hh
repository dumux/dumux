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
#ifndef DUMUX_RANS1P_SUBPROBLEM_HH
#define DUMUX_RANS1P_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

#include <dumux/freeflow/rans/twoeq/komega/model.hh>
#include <dumux/freeflow/rans/twoeq/komega/problem.hh>

namespace Dumux {

template <class TypeTag>
class FreeFlowSubProblem;

namespace Properties {

NEW_TYPE_TAG(RANSTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, KOmega));

// the fluid system
SET_PROP(RANSTypeTag, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::OnePGas<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
SET_TYPE_PROP(RANSTypeTag, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(RANSTypeTag, Problem, Dumux::FreeFlowSubProblem<TypeTag> );

SET_BOOL_PROP(RANSTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(RANSTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(RANSTypeTag, EnableGridVolumeVariablesCache, true);
}

/*!
 * \brief The free-flow sub problem
 */
template <class TypeTag>
class FreeFlowSubProblem : public KOmegaProblem<TypeTag>
{
    using ParentType = KOmegaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

public:
    FreeFlowSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "RANS"), eps_(1e-6), couplingManager_(couplingManager)
    {
        inletReynoldsNumber_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletReynoldsNumber");
        inletPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletPressure");
        temperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Temperature");

        darcyStart_ = getParamFromGroup<Scalar>(this->paramGroup(), "Grid.DarcyStart");
        darcyEnd_ = getParamFromGroup<Scalar>(this->paramGroup(), "Grid.DarcyEnd");
        smoothingZoneDistance_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.SmoothingZoneDistance");
        smoothingSlipVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.SmoothingSlipVelocity");

        Dumux::TurbulenceProperties<Scalar, FVGridGeometry::GridView::dimensionworld, true> turbulenceProperties;
        FluidState fluidState;
        fluidState.setPressure(0, inletPressure_);
        fluidState.setTemperature(temperature_);
        Scalar density = FluidSystem::density(fluidState, 0);
        Scalar kinematicViscosity = FluidSystem::viscosity(fluidState, 0) / density;
        Scalar charLength = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        inletVelocity_ = inletReynoldsNumber_ * kinematicViscosity / charLength;
        turbulentKineticEnergy_ = getParam<Scalar>("Problem.InletTurbulentKineticEnergy",
                                                   turbulenceProperties.turbulentKineticEnergy(inletVelocity_, charLength, kinematicViscosity, true));
        dissipation_ = getParam<Scalar>("Problem.InletDissipationRate",
                                        turbulenceProperties.dissipationRate(inletVelocity_, charLength, kinematicViscosity, true));
        std::cout << std::endl;
    }

    /*!
     * \name Boundary Locations
     */
    // \{

    bool isOnWall(const SubControlVolumeFace& scvf) const
    {
        GlobalPosition globalPos = scvf.ipGlobal();
        return isOnWallAtPos(globalPos);
    }

    bool isOnCouplingWall(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf);
    }

    bool isOnWallAtPos(const GlobalPosition& globalPos) const
    {
        return (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }
    // \}

    bool isSmoothingZoneAtPos(const GlobalPosition& globalPos) const
    {
        return (onLowerBoundary_(globalPos) &&
               (((globalPos[0] > ((darcyStart_ - smoothingZoneDistance_) - eps_)) && (globalPos[0] < (darcyStart_ + eps_))) ||
               ((globalPos[0] < ((darcyEnd_ + smoothingZoneDistance_) + eps_)) && (globalPos[0] > (darcyEnd_ - eps_)))));
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return temperature_; }

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

        bTypes.setDirichlet(Indices::turbulentKineticEnergyIdx);
        bTypes.setDirichlet(Indices::dissipationIdx);

        if (onRightBoundary_(globalPos))
        {
            bTypes.setDirichlet(Indices::pressureIdx);
            bTypes.setOutflow(Indices::turbulentKineticEnergyEqIdx);
            bTypes.setOutflow(Indices::dissipationEqIdx);
        }
        else
        {
            bTypes.setDirichlet(Indices::velocityXIdx);
            bTypes.setDirichlet(Indices::velocityYIdx);
        }

        if (isOnCouplingWall(scvf))
        {
            bTypes.setCouplingNeumann(Indices::conti0EqIdx);
            bTypes.setCouplingNeumann(scvf.directionIndex());
            bTypes.setBJS(1 - scvf.directionIndex());
        }

        // set a fixed dissipation (omega) in one cell
        if (isOnWall(scvf))
            bTypes.setDirichletCell(Indices::dissipationIdx);

        return bTypes;
    }

     /*!
      * \brief Evaluate the boundary conditions for a dirichlet values at the boundary.
      *
      * \param element The finite element
      * \param scvf the sub control volume face
      * \note used for cell-centered discretization schemes
      */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        return initialAtPos(scvf.ipGlobal());
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
        using std::pow;
        unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto wallDistance = ParentType::wallDistance_[elementIdx];
        values[Indices::dissipationEqIdx] = 6.0 * ParentType::kinematicViscosity_[elementIdx]
                                            / (ParentType::betaOmega() * pow(wallDistance, 2));
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
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);
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
        PrimaryVariables values(0.0);

        values[Indices::pressureIdx] = inletPressure();
        values[Indices::velocityXIdx] = inletVelocity();

        if(onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos))
            values[Indices::velocityXIdx] = 0.0;

        values[Indices::turbulentKineticEnergyIdx] = turbulentKineticEnergy_;
        values[Indices::dissipationIdx] = dissipation_;

        if (isSmoothingZoneAtPos(globalPos))
        {
        values[Indices::velocityXIdx] = smoothingSlipVelocity_;
        }

        if (isOnWallAtPos(globalPos))
        {
            values[Indices::turbulentKineticEnergyIdx] = 0.0;
        }

        return values;
    }

    //! \brief Returns the inlet velocity.
    const Scalar inletVelocity() const
    { return inletVelocity_ ;}

    //! \brief Returns the inlet pressure.
    const Scalar inletPressure() const
    { return inletPressure_; }

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

    // the height of the free-flow domain
    const Scalar height_() const
    { return this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]; }

    Scalar eps_;

    Scalar inletReynoldsNumber_;
    Scalar inletVelocity_;
    Scalar inletPressure_;
    Scalar temperature_;
    Scalar turbulentKineticEnergy_;
    Scalar dissipation_;
    Scalar darcyEnd_;
    Scalar darcyStart_;
    Scalar smoothingZoneDistance_;
    Scalar smoothingSlipVelocity_;

    TimeLoopPtr timeLoop_;

    std::shared_ptr<CouplingManager> couplingManager_;

};
} //end namespace Dumux

#endif // DUMUX_RANS1P_SUBPROBLEM_HH
