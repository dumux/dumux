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
        useFixedPressureAndSaturationBoundary_ = getParam<bool>("Problem.UseFixedPressureAndSaturationBoundary", false);
        distributeByVolume_ = getParam<bool>("Problem.DistributeByVolume", true);
        pc_ = getParam<Scalar>("Problem.CapillaryPressure");
        nonWettingMassFlux_ = getParam<Scalar>("Problem.NonWettingMassFlux", 5e-8);
        logfile_.open("logfile_" + this->name() + ".txt");
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
        if (isOutletPore_(scv))
            bcTypes.setAllDirichlet();
        return bcTypes;
    }

    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1.0e5;
        values[snIdx] = 0.0;
        return values;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! Evaluate the source term for all phases within a given sub-control-volume.
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);

        // We fix the mass flux of non-wetting injection at inlet of pore-network
        // The total inlet mass flux is distributed according to the ratio
        // of pore volume at inlet to the total volumess
        if (!useFixedPressureAndSaturationBoundary_ && isInletPore_(scv))
        {
            // values[pwIdx] = 1e8*(elemVolVars[scv].pressure(0) - 1e5);
            if (distributeByVolume_)
                values[snIdx] = nonWettingMassFlux_ * ( scv.volume()/sumInletPoresVolume_ );
            else
                values[snIdx] = nonWettingMassFlux_ * std::pow((this->gridGeometry().poreInscribedRadius(scv.dofIndex())), 4.0)/sumInletPoresVolume_;
        }
        return values / scv.volume();
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

    //! Loop over the scv in the domain to calculate the sum volume of inner inlet pores
    void calculateSumInletVolume()
    {
        sumInletPoresVolume_ = 0.0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (isInletPore_(scv))
                {
                    if (distributeByVolume_)
                        sumInletPoresVolume_ += this->gridGeometry().poreVolume(scv.dofIndex())/this->gridGeometry().coordinationNumber(scv.dofIndex());
                    else
                        sumInletPoresVolume_ += std::pow(this->gridGeometry().poreInscribedRadius(scv.dofIndex()), 4.0);
                }
            }
        }
    }

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
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
            const auto invaded = fluxVarsCache.invaded();
            if constexpr (useThetaRegularization)
            {
                if(!invaded)
                {
                    using std::max;
                    auto pcEntry = this->spatialParams().pcEntry(element, elemVolVars);
                    auto upwindVolVars = insideVolVars;
                    if (insideVolVars.pressure(nPhaseIdx) < outsideVolVars.pressure(nPhaseIdx))
                        upwindVolVars = outsideVolVars;
                    auto dp = upwindVolVars.capillaryPressure() / pcEntry - 1.0;
                    // Use a regularization function for theta
                    return applyInternalThetaConstraint_(fvGeometry, reg_.eval(dp,invaded));
                }
                else
                {
                    auto pcSnapoff = this->spatialParams().pcSnapoff(element, elemVolVars);
                    auto upwindVolVars = insideVolVars;
                    if (insideVolVars.pressure(wPhaseIdx) < outsideVolVars.pressure(wPhaseIdx))
                        upwindVolVars = outsideVolVars;
                    using std::min; using std::abs;
                    auto dp = upwindVolVars.capillaryPressure() / abs(pcSnapoff) - sign(pcSnapoff);
                    // Use a regularization function for theta
                    return applyInternalThetaConstraint_(fvGeometry, reg_.eval(dp,invaded));
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
    bool useFixedPressureAndSaturationBoundary_;
    bool distributeByVolume_;
    Scalar pc_;
    Scalar nonWettingMassFlux_;
    Scalar sumInletPoresVolume_;
    std::ofstream logfile_;
    std::shared_ptr<CouplingManager> couplingManager_;
    Dumux::PoreNetwork::Throat::Regularization<Scalar> reg_;
};
} //end namespace Dumux

#endif
