// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Free-flow sub-problem for the coupled 1p2c_1p2c free-flow/pore-network-model test
 */

#ifndef DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_ONEPTWOC_PROBLEM_FREEFLOW_HH
#define DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_ONEPTWOC_PROBLEM_FREEFLOW_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1pnc/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief  Free-flow sub-problem for the coupled 1p2c_1p2c free-flow/pore-network-model test
 *         A two-dimensional Stokes flow region coupled to a pore-network model.
 */
template <class TypeTag, class BaseProblem>
class FreeFlowOnePNCTestProblem :  public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr bool useMoles = ModelTraits::useMoles();

public:
    FreeFlowOnePNCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, "FreeFlow")
    , couplingManager_(couplingManager)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        initialPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialPressure", 1e5);
        initialMoleFraction_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialMoleFraction", 0.001);
        onlyDiffusion_ = getParam<bool>("Problem.OnlyDiffusion", false);
        pressureDifference_ = onlyDiffusion_ ? 0.0 : getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PressureDifference", 10);

#if !ISOTHERMAL
        initialTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialTemperature", 273.15 + 20.0);
#endif
    }

    /*!
     * \name Problem parameters
     */
    // \{

    const std::string& name() const
    { return problemName_; }

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
        const auto& globalPos = scvf.center(); //avoid ambiguities at corners

        if constexpr (ParentType::isMomentumProblem())
        {
            if (couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::poreNetworkIndex, scvf))
            {
                values.setCouplingNeumann(Indices::momentumYBalanceIdx);
                values.setCouplingNeumann(Indices::momentumXBalanceIdx);
            }
            else if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
                values.setAllNeumann();
            else
                values.setAllDirichlet();
        }
        else
        {
            if (couplingManager_->isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex, scvf))
                values.setAllCouplingNeumann();
            else if (!onlyDiffusion_ && onLeftBoundary_(globalPos))
                values.setAllDirichlet();
            else
                values.setAllNeumann();
        }
        return values;
    }

    /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return initialAtPos(globalPos);
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
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        using SlipVelocityPolicy = NavierStokesSlipVelocity<typename GridGeometry::DiscretizationMethod,  NavierStokes::SlipConditions::BJ>;
        using FluxHelper = NavierStokesMomentumBoundaryFlux<typename GridGeometry::DiscretizationMethod, SlipVelocityPolicy>;

        if constexpr (ParentType::isMomentumProblem())
        {

            if (couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::poreNetworkIndex, scvf))
            {
                values[scvf.normalAxis()] += couplingManager_->momentumCouplingCondition(
                    CouplingManager::freeFlowMomentumIndex, CouplingManager::poreNetworkIndex,
                    fvGeometry, scvf, elemVolVars
                );

                values += FluxHelper::slipVelocityMomentumFlux(
                    *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache
                );
            }
            else if (onLeftBoundary_(globalPos))
            {
                values = FluxHelper::fixedPressureMomentumFlux(
                    *this, fvGeometry, scvf, elemVolVars,
                    elemFluxVarsCache, initialPressure_ + pressureDifference_, true /*zeroNormalVelocityGradient*/
                );
            }
            else if (onRightBoundary_(globalPos))
            {
                values = FluxHelper::fixedPressureMomentumFlux(
                    *this, fvGeometry, scvf, elemVolVars,
                    elemFluxVarsCache, initialPressure_, true /*zeroNormalVelocityGradient*/
                );
            }
        }
        else
        {
            if (couplingManager_->isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex, scvf))
            {
                const auto& massCouplingCondition = couplingManager_->massCouplingCondition(
                    CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex,
                    fvGeometry, scvf, elemVolVars
                );
                values[Indices::conti0EqIdx] = massCouplingCondition[Indices::conti0EqIdx];
                values[Indices::conti0EqIdx + 1] = massCouplingCondition[Indices::conti0EqIdx + 1];

#if !ISOTHERMAL
                values[Indices::energyEqIdx] = couplingManager_->energyCouplingCondition(CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex,
                    fvGeometry, scvf, elemVolVars);
#endif
            }
            else if (!onlyDiffusion_ && onRightBoundary_(globalPos))
            {
                using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
                DirichletValues outsideBoundaryPriVars = initialAtPos(globalPos);
                values = FluxHelper::scalarOutflowFlux(
                    *this, element, fvGeometry, scvf, elemVolVars, std::move(outsideBoundaryPriVars)
                );
            }
        }

        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume
     */
    template<class ElementVolumeVariables>
    Sources source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        auto source = Sources(0.0);

        return source;
    }

    // The following function defines the initial conditions
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values(0.0); //velocity is 0.0

        if constexpr (!ParentType::isMomentumProblem())
        {
            values[Indices::pressureIdx] = initialPressure_;
            values[Indices::conti0EqIdx +1] = initialMoleFraction_;
#if !ISOTHERMAL
            values[Indices::temperatureIdx] = initialTemperature_;
#endif
        }

        return values;
    }

    /*!
     * \brief Returns true if the scvf lies on a porous slip boundary
     */
    bool onSlipBoundary(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        assert(scvf.isFrontal());
        return scvf.boundary() && couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::poreNetworkIndex, scvf);
    }

    /*!
     * \brief Returns the beta value
     */
    Scalar betaBJ(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf, const GlobalPosition& tangentialVector) const
    {
        const Scalar radius = couplingManager_->coupledPoreInscribedRadius(fvGeometry, scvf);
        return 5.73 / radius; // this value only an approximation of wall friction is considered
    }

    /*!
     * \brief Returns the velocity in the porous medium (which is 0 by default according to Saffmann).
     */
    VelocityVector porousMediumVelocity(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        return couplingManager_->interfaceThroatVelocity(fvGeometry, scvf);
    }

    // \}

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    std::string problemName_;
    static constexpr Scalar eps_ = 1e-6;
    Scalar deltaP_;
    Scalar initialPressure_;
    Scalar pressureDifference_;
    Scalar initialMoleFraction_;
    bool onlyDiffusion_;
#if !ISOTHERMAL
    Scalar initialTemperature_;
#endif
    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif
