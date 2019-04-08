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
 * \ingroup BoundaryTests
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

namespace Dumux {
template <class TypeTag>
class StokesSubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct StokesOneP { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOneP> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOneP> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Vertical flow from top to bottom with a parabolic velocity profile.
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFaceVariables = typename GetPropType<TypeTag, Properties::GridFaceVariables>::LocalView;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Stokes"), eps_(1e-6), couplingManager_(couplingManager)
    {
        inletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Velocity");
        pressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
        temperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Temperature");
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
     * \brief Returns the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return temperature_; }

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
        GlobalPosition globalPos = scvf.center();

        // vertical flow (1p_2p)
//        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
//        {
//            values.setCouplingNeumann(Indices::conti0EqIdx);
//            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
//            values.setBJS(Indices::momentumXBalanceIdx);
//        }
//        else
//        {
//            values.setDirichlet(Indices::velocityXIdx);
//            values.setDirichlet(Indices::velocityYIdx);
//        }

        // horizontal (1p2c_2p2c)
        if (onUpperBoundary_(globalPos) || onLeftBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
        }

        if (onRightBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
        }

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBJS(Indices::momentumXBalanceIdx);
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
     *
     * For this method, the \a values variable stores primary variables.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = scvf.dofPosition();
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        FluidState fluidState;
        updateFluidStateForBC_(fluidState, elemVolVars[scv].pressure());

        const Scalar density = fluidState.density(0);

        if(onLeftBoundary_(globalPos))
        {
            values[Indices::conti0EqIdx] = -xVelocity_(globalPos) * density;
        }

        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
        }
        return values;
    }

    // \}

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
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        // vertical flow (1p_2p)
//        values[Indices::pressureIdx] = pressure_;
//        values[Indices::velocityYIdx] = 4.0 * inletVelocity_ * globalPos[0] * (this->fvGridGeometry().bBoxMax()[0] - globalPos[0])
//                                      / (this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0])
//                                      / (this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0]);
        // horizontal flow (1p2c_2p2c)
        FluidState fluidState;
        updateFluidStateForBC_(fluidState, pressure_);

        const Scalar density = FluidSystem::density(fluidState, 0);

        values[Indices::pressureIdx] = pressure_ + density*this->gravity()[1]*(globalPos[1] - this->fvGridGeometry().bBoxMin()[1]);
        values[Indices::velocityXIdx] = xVelocity_(globalPos);
        return values;
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
              for the Beavers-Joseph-Saffman boundary condition.
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    { return couplingManager().couplingData().interfacePermeability(element, scvf); }

    /*!
     * \brief Returns the alpha value required as input parameter for the
              Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    { return 1.0; }

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

    //! Updates the fluid state to obtain required quantities for IC/BC
    void updateFluidStateForBC_(FluidState& fluidState, const Scalar pressure) const
    {
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(0, pressure);
        fluidState.setSaturation(0, 1.0);

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
        const Scalar vmax = inletVelocity_;
        return  4 * vmax * (globalPos[1] - this->fvGridGeometry().bBoxMin()[1]) * (this->fvGridGeometry().bBoxMax()[1] - globalPos[1])
                / (height_() * height_());
    }

    // the height of the free-flow domain
    const Scalar height_() const
    { return this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]; }

    Scalar eps_;
    Scalar inletVelocity_;
    Scalar pressure_;
    Scalar temperature_;
    std::string problemName_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif // DUMUX_STOKES_SUBPROBLEM_HH
