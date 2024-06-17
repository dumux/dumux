// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for the two-phase pore network model.
 */
#ifndef DUMUX_PNM2P_PROBLEM_HH
#define DUMUX_PNM2P_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include "../regularization.hh"

namespace Dumux {

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        pwIdx = Indices::pressureIdx,
        snIdx = Indices::saturationIdx,
        wPhaseIdx = FluidSystem::phase0Idx,
        nPhaseIdx = FluidSystem::phase1Idx,
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    using GridFluxVariablesCache = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using InvasionState = std::decay_t<decltype(std::declval<GridFluxVariablesCache>().invasionState())>;
    static constexpr bool useThetaRegularization = InvasionState::stateMethod == Dumux::PoreNetwork::StateSwitchMethod::theta;

public:
    template<class SpatialParams>
    DrainageProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<SpatialParams> spatialParams,
                    const std::string& paramGroup = "",
                    std::shared_ptr<CouplingManager> couplingManager = nullptr)
    : ParentType(gridGeometry, spatialParams, paramGroup),
      couplingManager_(couplingManager),
      reg_(getParamFromGroup<std::string>(paramGroup, "Problem.RegularizationFunction", "Sine"),
           getParamFromGroup<Scalar>(paramGroup, "Problem.RegularizationDelta", 1e-8))
    {
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        pc_ = getParam<Scalar>("Problem.CapillaryPressure");
        logfile_.open("logfile_" + this->name() + ".txt");
        boundConductiviy_ = getParam<bool>("Problem.BoundConductivity", true);
        useUpwindPc_ = getParam<bool>("Problem.UseUpwindPc", true);
        switchConductivityCurve_ = getParam<bool>("Problem.switchConductivityCurve", true);
    }

     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        if (isInletPore_(scv) || isOutletPore_(scv))
           bcTypes.setAllDirichlet();
        else
           bcTypes.setAllNeumann();
        return bcTypes;
    }

    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        values[snIdx] = 0.0;
        // A global phase pressure difference (pn,inlet - pw,outlet) is specified and the saturation shall also be fixed, apply:
        // pw,inlet = pw,outlet = 1e5; pn,outlet = pw,outlet + pc(S=0) = pw,outlet; pn,inlet = pw,inlet + pc_
        if (isInletPore_(scv))
            values[snIdx] = 1.0 - this->spatialParams().fluidMatrixInteraction(element, scv, int()/*dummyElemsol*/).sw(pc_);
        return values;
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        values[snIdx] = 0.0;
        return values;
    }

    //!  Evaluate the initial invasion state of a pore throat
    bool initialInvasionState(const Element& element) const
    { return false; }
    // \}

    // throat parameter indicating invasion state
    template<class FluxVariablesCache, class SubControlVolumeFace>
    Scalar theta(const Element& element,
                 const FVElementGeometry& fvGeometry,
                 const ElementVolumeVariables& elemVolVars,
                 const FluxVariablesCache& fluxVarsCache,
                 const SubControlVolumeFace& scvf) const
    {
        if constexpr (std::is_void_v<CouplingManager>)
        {
            auto upwindW = (elemVolVars[0].pressure(wPhaseIdx) > elemVolVars[1].pressure(wPhaseIdx)) ?
                            elemVolVars[0] : elemVolVars[1];
            auto upwindN = (elemVolVars[0].pressure(nPhaseIdx) > elemVolVars[1].pressure(nPhaseIdx)) ?
                            elemVolVars[0] : elemVolVars[1];
            const auto invaded = fluxVarsCache.invaded();
            if constexpr (useThetaRegularization)
            {
                if (switchConductivityCurve_)
                {
                    if(!invaded)
                    {
                        using std::max;
                        auto pcEntry = this->spatialParams().pcEntry(element, elemVolVars);
                        auto dp = max(elemVolVars[0].capillaryPressure(),
                                      elemVolVars[1].capillaryPressure()) / pcEntry - 1.0;
                        if (useUpwindPc_)
                            dp = upwindN.capillaryPressure()/pcEntry - 1.0;
                        // Use a regularization function for theta
                        return applyInternalThetaConstraint_(fvGeometry, reg_.eval(dp,invaded));
                    }
                    else
                    {
                        using std::min; using std::abs;
                        auto pcSnapoff = this->spatialParams().pcSnapoff(element, elemVolVars);
                        auto dp = min(elemVolVars[0].capillaryPressure(),
                                      elemVolVars[1].capillaryPressure()) / abs(pcSnapoff) - sign(pcSnapoff);
                        if (useUpwindPc_)
                            dp = upwindW.capillaryPressure() / abs(pcSnapoff) - sign(pcSnapoff);
                        if (boundConductiviy_ && upwindN.saturation(nPhaseIdx) < 1e-3)
                            return M_PI + upwindN.saturation(nPhaseIdx);
                        // Use a regularization function for theta
                        return applyInternalThetaConstraint_(fvGeometry, reg_.eval(dp,invaded));
                    }
                }
                else
                {
                    using std::max;
                    auto pcEntry = this->spatialParams().pcEntry(element, elemVolVars);
                    auto dp = max(elemVolVars[0].capillaryPressure(),
                                  elemVolVars[1].capillaryPressure()) / pcEntry - 1.0;
                    if (useUpwindPc_)
                        dp = upwindN.capillaryPressure()/pcEntry - 1.0;
                    // Use a regularization function for theta
                    return applyInternalThetaConstraint_(fvGeometry, reg_.eval(dp, false));
                }
            }
            else
            {
                Scalar val = invaded ? 1.0 : 0.0;
                return applyInternalThetaConstraint_(fvGeometry, val);
            }
        }
        else
            return couplingManager_->theta(element);
    }

    // define when vtk output is written
    bool shouldWriteOutput(const int timeStepIndex, const GridVariables& gridVariables) const
    {
        if (vtpOutputFrequency_ < 0)
            return true;

        if (vtpOutputFrequency_ == 0)
            return (timeStepIndex == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
        else
            return (timeStepIndex % vtpOutputFrequency_ == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
    }

    /*!
     * \brief Called at the end of each time step
     */
    template<class AveragedValues>
    void postTimeStep(const Scalar time, const AveragedValues& avgValues, std::size_t numThroatsInvaded, const Scalar dt)
    {
        logfile_ << std::fixed << std::left << std::setw(20) << std::setfill(' ') << time
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgSat"]
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPw"]
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"]
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"] - avgValues["avgPw"]
                 << std::left << std::setw(20) << std::setfill(' ') << numThroatsInvaded
                 << std::endl;
    }

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

    Scalar applyInternalThetaConstraint_(const FVElementGeometry& fvGeometry,
                                         Scalar thetaUnconstrained) const
    {
        static const auto blockNonwettingPhase = getParamFromGroup<std::vector<int>>(this->paramGroup(), "InvasionState.BlockNonwettingPhaseAtThroatLabel", std::vector<int>{});
        auto theta = thetaUnconstrained;

        if(!blockNonwettingPhase.empty())
        {
            auto eIdx = fvGeometry.gridGeometry().elementMapper().index(fvGeometry.element());
            if(std::find(blockNonwettingPhase.begin(), blockNonwettingPhase.end(), fvGeometry.gridGeometry().throatLabel(eIdx)) != blockNonwettingPhase.end())
                theta = 0.0;
        }

        return theta;
    }

    int vtpOutputFrequency_;
    Scalar pc_;
    std::ofstream logfile_;
    std::shared_ptr<CouplingManager> couplingManager_;
    Dumux::PoreNetwork::Throat::Regularization<Scalar> reg_;
    bool useUpwindPc_;
    bool boundConductiviy_;
    bool switchConductivityCurve_;
};
} //end namespace Dumux

#endif
