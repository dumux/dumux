// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Free-flow sub-problem for the coupled 1p_1p free-flow/pore-network-model test
 */

#ifndef DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_PROBLEM_FREEFLOW_HH
#define DUMUX_TEST_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORE_NETWORK_PROBLEM_FREEFLOW_HH

#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief  Free-flow sub-problem for the coupled 1p_1p free-flow/pore-network-model test
 *         A two-dimensional Stokes flow region coupled to a pore-network model.
 */
template <class TypeTag, class BaseProblem>
class FreeFlowOnePTestProblem :  public BaseProblem
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

public:
    FreeFlowOnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, "FreeFlow")
    , couplingManager_(couplingManager)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        deltaP_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PressureDifference", 0.0);
        outletPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.OutletPressure", 1e5);
        deltaT_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.TemperatureDifference", 10.0);
        outletTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.OutletTemperature", 273.15 + 20.0);
        verticalFlow_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.VerticalFlow", false);
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
                values.setCouplingNeumann(Indices::momentumXBalanceIdx);
                values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            }
            else if (!verticalFlow_ && (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)))
                values.setAllNeumann(); //inflow/outflow
            else if (verticalFlow_ && onUpperBoundary_(globalPos))
                values.setAllNeumann(); //outflow
            else
                values.setAllDirichlet(); //e.g. velocities at walls (parallel to flow direction)
        }
        else
        {
            if (couplingManager_->isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex, scvf))
                values.setAllCouplingNeumann(); //mass and energy coupling
            else
            {
                if (!verticalFlow_ && (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)))
                {
                    values.setDirichlet(Indices::pressureIdx);
                    values.setDirichlet(Indices::temperatureIdx);
                }
                else
                    values.setAllNeumann();
            }
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
        DirichletValues values = initialAtPos(globalPos);

        if constexpr (!ParentType::isMomentumProblem())
        {
            if (onLeftBoundary_(globalPos))
            {
                values[Indices::pressureIdx] = outletPressure_ + deltaP_;
                values[Indices::temperatureIdx] = outletTemperature_ + deltaT_;
            }
        }
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment/ control volume.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
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
        using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;

        if constexpr (ParentType::isMomentumProblem())
        {

            if (couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::poreNetworkIndex, scvf))
            {
                values += couplingManager_->momentumCouplingCondition(
                    CouplingManager::freeFlowMomentumIndex, CouplingManager::poreNetworkIndex,
                    fvGeometry, scvf, elemVolVars
                );

                values += FluxHelper::slipVelocityMomentumFlux(
                    *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache
                );
            }
            else //inlet, outlet
            {
                const Scalar inletPressure = onLeftBoundary_(globalPos) ? outletPressure_ + deltaP_ : outletPressure_;
                values = FluxHelper::fixedPressureMomentumFlux(
                    *this, fvGeometry, scvf, elemVolVars,
                    elemFluxVarsCache, inletPressure, true /*zeroNormalVelocityGradient*/
                );
            }
        }
        else
        {
            if (couplingManager_->isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex, scvf))
            {
                values[Indices::conti0EqIdx] = couplingManager_->massCouplingCondition(
                    CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex,
                    fvGeometry, scvf, elemVolVars);
                values[Indices::energyEqIdx] = couplingManager_->energyCouplingCondition(CouplingManager::freeFlowMassIndex, CouplingManager::poreNetworkIndex,
                    fvGeometry, scvf, elemVolVars);
            }
            else
            {
                using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
                DirichletValues outsideBoundaryPriVars;
                outsideBoundaryPriVars[Indices::pressureIdx] = outletPressure_;
                outsideBoundaryPriVars[Indices::temperatureIdx] = outletTemperature_;
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
            values[Indices::pressureIdx] = outletPressure_;
            values[Indices::temperatureIdx] = outletTemperature_;
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
        return 5.73 / radius; // this values is only an approximation of wall friction is considered
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
    Scalar outletPressure_;
    Scalar deltaT_;
    Scalar outletTemperature_;
    bool verticalFlow_;

    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif