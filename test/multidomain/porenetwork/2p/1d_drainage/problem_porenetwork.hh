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
#ifndef DUMUX_PNM2P_MD_PROBLEM_HH
#define DUMUX_PNM2P_MD_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

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
        nPhaseIdx = FluidSystem::phase1Idx,

#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif
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
      couplingManager_(couplingManager)
    {
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        useFixedPressureAndSaturationBoundary_ = getParam<bool>("Problem.UseFixedPressureAndSaturationBoundary", false);
        distributeByVolume_ = getParam<bool>("Problem.DistributeByVolume", true);
        pc_ = getParam<Scalar>("Problem.CapillaryPressure");
        nonWettingMassFlux_ = getParam<Scalar>("Problem.NonWettingMassFlux", 5e-8);
        logfile_.open("time_steps_" + this->name() + ".txt");
        regDelta_ = getParamFromGroup<Scalar>(paramGroup, "Problem.RegularizationDelta", 1e-8);
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \name Problem parameters
     */
    // \{

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
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
#if DRAINAGE
        if (isOutletPore_(scv))
            bcTypes.setAllDirichlet();
#else
        bcTypes.setAllDirichlet();
#endif
        return bcTypes;
    }


    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
#if DRAINAGE
        values[pwIdx] = 1.0e5;
        values[snIdx] = 0.0;
#else
        if (isOutletPore_(scv))
        {
            values[pwIdx] = 1.0e5;
            values[snIdx] = 0.0;
        }
        else if (isInletPore_(scv))
        {
            values[pwIdx] = 1.0e5;
            values[snIdx] = 1.0;
        }
#endif

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
        // The total inlet mass flux is distributed according to the ratio of pore volume
        if (isInletPore_(scv))
            values[snIdx] = nonWettingMassFlux_;
        return values / scv.volume();
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        values[snIdx] = 0.0;
#if !DRAINAGE
        const auto dofIdxGlobal = this->gridGeometry().vertexMapper().index(vertex);
        if (isInletPore_(dofIdxGlobal))
        {
            values[pwIdx] = 1e5;
            values[snIdx] = 1.0;
        }
#endif
        return values;
    }



    //!  Evaluate the initial invasion state of a pore throat
#if DRAINAGE
    bool initialInvasionState(const Element& element) const
    { return false; }
#else
    bool initialInvasionState(const Element& element) const
    { return true; }
#endif
    // \}


    template<class FluxVariablesCache>
    Scalar theta(const Element& element,
                 const FVElementGeometry& fvGeometry,
                 const ElementVolumeVariables& elemVolVars,
                 const FluxVariablesCache& fluxVarsCache) const
    {
        if constexpr (std::is_void_v<CouplingManager>)
        {
            const auto invaded = fluxVarsCache.invaded();
            if constexpr (useThetaRegularization)
            {
                if(!invaded)
                {
                    using std::max;
                    auto pcEntry = this->spatialParams().pcEntry(element, elemVolVars);
                    auto dp = max(elemVolVars[0].capillaryPressure(),
                                  elemVolVars[1].capillaryPressure()) / pcEntry - 1.0;
                    // Use a regularized heavyside function for theta
                    auto theta = regHeaviside_(dp,invaded);

                    return std::min(std::max(0.0,theta),1.0);
                }
                else
                {
                    using std::min; using std::abs;
                    auto pcSnapoff = this->spatialParams().pcSnapoff(element, elemVolVars);
                    auto dp = min(elemVolVars[0].capillaryPressure(),
                                  elemVolVars[1].capillaryPressure()) / abs(pcSnapoff) - sign(pcSnapoff);
                    // Use a regularized heavyside function for theta
                    auto theta = regHeaviside_(dp, invaded);

                    return std::min(std::max(0.0,theta),1.0);
                }
            }
            else
            {
                return invaded ? 1.0 : 0.0;
            }
        }
        else
            return couplingManager_->theta(element);
    }

    // \}

    /*!
     * \brief Called at the end of each time step
     */
    void writeInvasionTimeStep(const Scalar time)
    {
        logfile_ << std::fixed << std::left << std::setw(20)
                 << std::setfill(' ') << time << std::endl;
    }

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == 2;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == 1;
    }

    Scalar regAbs_(Scalar v) const
    {
        return v*v/(std::sqrt(v*v + regDelta_*regDelta_));
    }

    Scalar regSign_(Scalar v) const
    {
        return v/(regAbs_(v) + regDelta_);
    }

    Scalar regHeaviside_(Scalar v, bool invaded) const
    {
        //  0.5*(1 + regSign_(dp))
        if(!invaded)
        {
            using std::sin; using std::min; using std::max;
            v = max(0.0,min(v,2*regDelta_));
            return 0.5*(1+sin(M_PI*(v-regDelta_)/(2.0*regDelta_)));
        }
        else
        {
            using std::sin; using std::min; using std::max;
            v = max(-2*regDelta_,min(v,0.0));
            return 0.5*(1+sin(M_PI*(v+regDelta_)/(2.0*regDelta_)));
        }
    }

    int vtpOutputFrequency_;
    bool useFixedPressureAndSaturationBoundary_;
    bool distributeByVolume_;
    Scalar pc_;
    Scalar nonWettingMassFlux_;
    Scalar sumInletPoresVolume_;
    Scalar regDelta_;
    std::ofstream logfile_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace Dumux

#endif