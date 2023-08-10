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

#include <memory>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetwork/common/outletpcgradient.hh>

namespace Dumux {

template <class TypeTag>
class DrainageProblem;

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
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
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

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
    using OutletCapPressureGradient = typename Dumux::PoreNetwork::OutletCapPressureGradient<GridVariables, SolutionVector>;

public:
    template<class SpatialParams>
    DrainageProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        initialPc_ = getParam<Scalar>("Problem.InitialPc", 3000);
        finalPc_ = getParam<Scalar>("Problem.FinalPc", 7000);
        numSteps_ = getParam<int>("Problem.NumSteps", 10);
        swShiftThreshold_ = getParam<Scalar>("Problem.RelShiftThreshold", 1e-6);
        writeOnlyEqPoints_ = getParam<bool>("Problem.WriteOnlyEquilibriumPoints", false);

        pcEpisopde_.resize(numSteps_ + 1);
        for (int i = 0 ; i < pcEpisopde_.size(); i++)
              pcEpisopde_[i] = initialPc_ + i*(finalPc_ - initialPc_)/numSteps_;

        std::cout << "The following global PCs are applied: " << std::endl;
        for (auto x: pcEpisopde_)
        {
            std::cout << x << std::endl;
        }

        if (!writeOnlyEqPoints_)
        {
            logfile_.open("logfile_" + this->name() + ".txt"); //for the logfile
            logfile_ <<"Logfile for: " + this->name()  << std::endl;
            logfile_ << std::left << std::setw(20) << std::setfill(' ') << "Time"
                     << std::left << std::setw(20) << std::setfill(' ') << "globalPc"
                     << std::left << std::setw(20) << std::setfill(' ') << "swAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "pwAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "pnAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "pcAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "numThroatsInvaded"
                     << std::endl;
        }

        logfileEqPoints_.open("eqPoints_" + this->name() + ".txt");
        logfileEqPoints_ << std::left << std::setw(20) << std::setfill(' ') << "Time"
                         << std::left << std::setw(20) << std::setfill(' ') << "globalPc"
                         << std::left << std::setw(20) << std::setfill(' ') << "swAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "pwAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "pnAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "pcAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "numThroatsInvaded"
                         << std::endl;
        step_ = 0;
    }


    /*!
     * \brief Called at the end of each time step
     */
    template<class AveragedValues>
    void postTimeStep(const Scalar time, const AveragedValues& avgValues, std::size_t numThroatsInvaded, const Scalar dt)
    {
        const Scalar avgSw = avgValues["avgSat"];

        if (!writeOnlyEqPoints_)
        {
            logfile_ << std::fixed << std::left << std::setw(20) << std::setfill(' ') << time
                                   << std::left << std::setw(20) << std::setfill(' ') <<  pcEpisopde_[step_]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgSat"]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPw"]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"] - avgValues["avgPw"]
                                   << std::left << std::setw(20) << std::setfill(' ') << numThroatsInvaded
                                   << std::endl;
        }

        // store the three most recent averaged saturations
        std::rotate(swAvg_.rbegin(), swAvg_.rbegin()+1, swAvg_.rend());
        swAvg_[0]= avgSw;

        // Check for steady state and end episode
        dSwDt_ = std::abs(swAvg_[0]-swAvg_[1])/dt;
        const Scalar pc = pcEpisopde_[step_];
        std::cout << "global pC applied: " << pc << " / " << finalPc_ << " (step " << step_ << " of " << numSteps_ << ")" << std::endl;
        std::cout << "swAverage: " << swAvg_[0] << " (relative shift: " << dSwDt_ << "). " << std::endl;
        std::cout << numThroatsInvaded << " of " << this->gridGeometry().gridView().size(0) << " throats invaded." << std::endl;
        if(dSwDt_ < swShiftThreshold_)
        {
            std::cout << "Equlibrium point reached!" << std::endl;

            logfileEqPoints_ << std::fixed << std::left << std::setw(20) << std::setfill(' ') << time
                                           << std::left << std::setw(20) << std::setfill(' ') <<  pcEpisopde_[step_]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgSat"]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPw"]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"] - avgValues["avgPw"]
                                           << std::left << std::setw(20) << std::setfill(' ') << numThroatsInvaded
                                           << std::endl;
            inEquilibrium_ = true;
            ++step_;
        }
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

    bool equilibriumPointReached() const
    { return dSwDt_ < swShiftThreshold_; }

     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume &scv) const
    {
        BoundaryTypes bcTypes;
        if (isInletPore_(scv) || isOutletPore_(scv))
           bcTypes.setAllDirichlet();

        else // neuman for the remaining boundaries
           bcTypes.setAllNeumann();

        return bcTypes;
    }


    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        values[snIdx] = 0.0;

        const auto& fluidMatrixInteraction = this->spatialParams().fluidMatrixInteraction(element, scv, 0);
        if (isInletPore_(scv))
            values[snIdx] = 1.0 - fluidMatrixInteraction.sw(pcEpisopde_[step_]);
        else if (isOutletPore_(scv))
        {
            values[snIdx] = 1.0 - outletPcGradient_->zeroPcGradientSw(element, scv);
        }

        return values;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! Evaluate the source term for all phases within a given sub-control-volume.
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        return PrimaryVariables(0.0);
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        return values;
    }

    //!  Evaluate the initial invasion state of a pore throat
    bool initialInvasionState(const Element& element) const
    { return false; }

    // \}

    void outletCapPressureGradient(std::shared_ptr<OutletCapPressureGradient> outletPcGradient)
    {  outletPcGradient_ = outletPcGradient;}

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

    int vtpOutputFrequency_;
    Scalar initialPc_;
    Scalar finalPc_;
    int numSteps_;
    std::vector<Scalar> pcEpisopde_;
    std::array<Scalar, 2> swAvg_ = {{1.0, 1.0}};
    std::ofstream logfile_;
    std::ofstream logfileEqPoints_;
    Scalar dSwDt_ = 0.0;
    mutable bool inEquilibrium_ = false;
    Scalar swShiftThreshold_;
    bool writeOnlyEqPoints_;
    std::shared_ptr<OutletCapPressureGradient> outletPcGradient_;

    int step_;
};
} //end namespace Dumux

#endif
