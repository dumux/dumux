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
 *
 * \brief A test problem for the two-phase pore network model.
 */
#ifndef DUMUX_PNM2P_PROBLEM_HH
#define DUMUX_PNM2P_PROBLEM_HH

#include <memory>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template <class TypeTag>
class DrainageProblem;

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        pwIdx = Indices::pressureIdx,
        snIdx = Indices::saturationIdx,
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
      couplingManager_(couplingManager)
    {
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        useFixedPressureAndSaturationBoundary_ = getParam<bool>("Problem.UseFixedPressureAndSaturationBoundary", false);
        distributeByVolume_ = getParam<bool>("Problem.DistributeByVolume", true);
        pc_ = getParam<Scalar>("Problem.CapillaryPressure");
        nonWettingMassFlux_ = getParam<Scalar>("Problem.NonWettingMassFlux", 5e-8);
        logfile_.open("logfile_" + this->name() + ".txt");
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

        // If a global phase pressure difference (pn,inlet - pw,outlet) with fixed saturations is specified, use a Dirichlet BC here
        if (useFixedPressureAndSaturationBoundary_ && isInletPore_(scv))
            bcTypes.setAllDirichlet();
        else if (isOutletPore_(scv))
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

        // If a global phase pressure difference (pn,inlet - pw,outlet) is specified and the saturation shall also be fixed, apply:
        // pw,inlet = pw,outlet = 1e5; pn,outlet = pw,outlet + pc(S=0) = pw,outlet; pn,inlet = pw,inlet + pc_
        if (useFixedPressureAndSaturationBoundary_ && isInletPore_(scv))
            values[snIdx] = 1.0 - this->spatialParams().fluidMatrixInteraction(element, scv, int()/*dummyElemsol*/).sw(pc_);
        else if (isOutletPore_(scv))
        {
            values[pwIdx] = 1.0e5;
            values[snIdx] = 0.0;
        }
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

    /*!
     * \brief Called at the end of each time step
     */
    template<class AveragedValues>
    void postTimeStep(const Scalar time, const AveragedValues& avgValues, std::size_t numThroatsInvaded, const Scalar dt)
    {
        const Scalar avgSw = avgValues["avgSat"];


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
    Scalar regAbs_(Scalar v) const
    {
        return v*v/(std::sqrt(v*v + regDelta_*regDelta_));
    }

    Scalar regSign_(Scalar v) const
    {
        return v/(regAbs_(v) + regDelta_);
    }

    Scalar regHeaviside_(Scalar v) const
    {
        //  0.5*(1 + regSign_(dp))
        using std::sin; using std::min; using std::max;
        v = max(0.0,min(v,2*regDelta_));
        return 0.5*(1+sin(M_PI*(v-regDelta_)/(2.0*regDelta_)));
    }

    int vtpOutputFrequency_;
    bool useFixedPressureAndSaturationBoundary_;
    bool distributeByVolume_;
    Scalar pc_;
    Scalar regDelta_;
    Scalar nonWettingMassFlux_;
    Scalar sumInletPoresVolume_;
    std::ofstream logfile_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace Dumux

#endif
