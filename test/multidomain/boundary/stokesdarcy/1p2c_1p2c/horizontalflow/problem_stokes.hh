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
 * \brief A simple Navier-Stokes test problem for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

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
// Create new type tags
namespace TTag {
struct StokesOnePTwoC { using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>; };
} // end namespace TTag

// The fluid system
SET_PROP(StokesOnePTwoC, FluidSystem)
{
  using H2OAir = FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>;
  static constexpr auto phaseIdx = H2OAir::liquidPhaseIdx; // simulate the water phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// Set the grid type
SET_TYPE_PROP(StokesOnePTwoC, Grid, Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(StokesOnePTwoC, Problem, Dumux::StokesSubProblem<TypeTag> );

SET_BOOL_PROP(StokesOnePTwoC, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(StokesOnePTwoC, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(StokesOnePTwoC, EnableGridVolumeVariablesCache, true);

// Use moles
SET_BOOL_PROP(StokesOnePTwoC, UseMoles, true);

// Do not replace one equation with a total mass balance
SET_INT_PROP(StokesOnePTwoC, ReplaceCompEqIdx, 3);
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
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Stokes"), eps_(1e-6), injectionState_(false), couplingManager_(couplingManager)
    {
        inletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Velocity");
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
        inletMoleFraction_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletMoleFraction");
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

    bool shouldWriteRestartFile() const
    { return false; }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

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

        const auto& globalPos = scvf.dofPosition();

        if(onLeftBoundary_(globalPos))
        {
            values.setDirichlet(Indices::conti0EqIdx + 1);
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else if(onRightBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(Indices::conti0EqIdx + 1);
        }
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(Indices::conti0EqIdx + 1);
        }

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
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
     * \param scvf The sub control volume face
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(globalPos);

        // start injecting after the velocity field had enough time to initialize
        if(globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_ && isInjectionPeriod())
            values[Indices::conti0EqIdx + 1] = inletMoleFraction_;

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
        NumEqVector values(0.0);

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto tmp = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::conti0EqIdx] = tmp[0];
            values[Indices::conti0EqIdx + 1] = tmp[1];
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
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = pressure_;
        values[Indices::velocityXIdx] = inletVelocity_ * (globalPos[1] - this->fvGridGeometry().bBoxMin()[1])
                                              * (this->fvGridGeometry().bBoxMax()[1] - globalPos[1])
                                              / (0.25 * (this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1])
                                              * (this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]));

        return values;
    }

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
        return 1.0;
    }

    void setInjectionState(const bool yesNo)
    {
        injectionState_ = yesNo;
    }

    bool isInjectionPeriod() const
    {
        return injectionState_;
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

    Scalar eps_;
    Scalar inletVelocity_;
    Scalar pressure_;
    Scalar inletMoleFraction_;
    bool injectionState_;
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace

#endif // DUMUX_STOKES_SUBPROBLEM_HH
