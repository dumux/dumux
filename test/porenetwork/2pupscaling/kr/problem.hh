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
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

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
        pressure_ = getParam<Scalar>("Problem.Pressure", 0.0);
        saturationw_ = getParam<Scalar>("Problem.SaturationWetting", 1.0);
        minSource_ = getParam<Scalar>("Problem.MinSource", 0.0);
        maxSource_ = getParam<Scalar>("Problem.MaxSource", 1e-12);

        useLabels_ = getParam<bool>("Problem.UseLabels", true);
        eps_ = getParam<Scalar>("Problem.Epsilon", 1e-7);

        sourceEpisopdeN_.resize(numSteps_ + 1);
        for (int i = 0 ; i < sourceEpisopdeN_.size(); i++)
              sourceEpisopdeN_[i] = minSource_ + i*(maxSource_ - minSource_) / numSteps_;

        sourceEpisopdeW_ = sourceEpisopdeN_;
        std::reverse(sourceEpisopdeW_.begin(), sourceEpisopdeW_.end());

        for (int i = 0 ; i < sourceEpisopdeN_.size(); i++)
            std::cout<<" sourceEpisopdeN_[i] "<<sourceEpisopdeN_[i]<<" sourceEpisopdeW_[i] "<<sourceEpisopdeW_[i]<<std::endl;

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
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume &scv) const
    {
        BoundaryTypes bcTypes;
        if (isInletPore_(scv))
            bcTypes.setAllNeumann();
        else if (isOutletPore_(scv))
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
        if (isOutletPore_(scv))
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
        PrimaryVariables values(0.0);
        if (isInletPore_(scv))
        {
            const auto& volVar =  elemVolVars[scv.indexInElement()];
            values[Indices::conti0EqIdx] = /*sourceEpisopdeW_[step_]*/0.0* volVar.density(FluidSystem::phase0Idx);
            values[Indices::conti0EqIdx + 1] = /*sourceEpisopdeN_[step_]*/1e-9* volVar.density(FluidSystem::phase1Idx);
        }
        values /= scv.volume();

        return values;
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = pressure_;
        values[snIdx] = 1.0 - saturationw_;
        return values;
    }

    //!  Evaluate the initial invasion state of a pore throat
    bool initialInvasionState(const Element& element) const
    {   if (saturationw_ == 0.0)
            return true;
        return false;
    }

    // \}

    void outletCapPressureGradient(std::shared_ptr<OutletCapPressureGradient> outletPcGradient)
    {  outletPcGradient_ = outletPcGradient;}

    // Return the label of inlet pores assuming a previously set direction.
    int inletPoreLabel() const
    {
        static constexpr std::array<int, 3> label = {2, 2, 2};//{1, 3, 5};
        return label[direction_];
    }

    // Return the label of outlet pores assuming a previously set direction.
    int outletPoreLabel() const
    {
        static constexpr std::array<int, 3> label = {3, 3, 3};//{2, 4, 6};
        return label[direction_];
    }

    // Set the current direction (0:x, 1:y, 2:z) in which the pressure gradient is applied
    void setDirection(int directionIdx)
    { direction_ = directionIdx; }

    // Get the current direction in which the pressure gradient is applied.
    int direction() const
    { return direction_; }

    void nextStep()
    {   step_++; }

    bool finished()
    {   return step_>numSteps_; }

private:

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == Labels::inlet;
    }


    bool isInletPore_(const SubControlVolume& scv) const
    {
        if (useLabels_)
            return inletPoreLabel() == this->gridGeometry().poreLabel(scv.dofIndex());
        else
            return scv.dofPosition()[direction_] < this->gridGeometry().bBoxMin()[direction_] + eps_;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        if (useLabels_)
            return outletPoreLabel() == this->gridGeometry().poreLabel(scv.dofIndex());
        else
            return scv.dofPosition()[direction_] > this->gridGeometry().bBoxMax()[direction_] - eps_;
    }

    int vtpOutputFrequency_;
    Scalar initialPc_;
    Scalar finalPc_;
    int numSteps_;
    std::ofstream logfile_;
    std::ofstream logfileEqPoints_;
    Scalar dSwDt_ = 0.0;
    mutable bool inEquilibrium_ = false;
    Scalar swShiftThreshold_;
    bool writeOnlyEqPoints_;
    std::shared_ptr<OutletCapPressureGradient> outletPcGradient_;
    Scalar pressure_;
    Scalar saturationw_;
    int direction_;
    bool useLabels_;
    Scalar eps_;
    Scalar minSource_;
    Scalar maxSource_;
    std::vector<Scalar> sourceEpisopdeW_;
    std::vector<Scalar> sourceEpisopdeN_;
    int step_ = 0;
};
} //end namespace Dumux

#endif
