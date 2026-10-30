// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Driver for the quasi-static two-phase invasion-percolation sweep, run once as a
 *        preprocessing step to obtain the per-throat, per-saturation-step phase conductances
 *        used by the single-phase upscaling solves.
 */
#ifndef DUMUX_PNM_TWOP_UPSCALING_PROBLEM_STATIC_HH
#define DUMUX_PNM_TWOP_UPSCALING_PROBLEM_STATIC_HH

#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <fstream>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/porenetwork/common/throatproperties.hh>
#include <dumux/flux/porenetwork/advection.hh> // for Detail::Transmissibility

#include <dumux/porenetwork/2p/static/staticpnm.hh>
#include <dumux/io/gnuplotinterface.hh>

namespace Dumux::PoreNetwork {

/*!
 * \brief Stand-in types for the Problem/ElementVolumeVariables arguments expected by the
 *        generic single-phase throat-transmissibility interface.
 *
 * `Detail::Transmissibility<...>::singlePhaseTransmissibility(problem, element, fvGeometry,
 * scvf, elemVolVars, fluxVarsCache, phaseIdx)` is only ever called here with
 * `TransmissibilityPatzekSilin<Scalar>` using its default `considerPoreResistance = false`.
 * In that case, the only code path that would actually dereference `problem`/`elemVolVars` is
 * compiled away by `if constexpr (!considerPoreResistance)` inside
 * `TransmissibilityPatzekSilin::singlePhaseTransmissibility` (see
 * dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh), so these
 * placeholders are never actually used and only need to exist as valid types. If
 * `considerPoreResistance` were ever switched on, real `Problem`/`ElementVolumeVariables`
 * objects would be required here instead.
 */
struct StaticUpscalingMockProblem {};
struct StaticUpscalingMockElementVolumeVariables {};

/*!
 * \ingroup PoreNetworkModels
 * \brief Runs the quasi-static two-phase invasion-percolation model (TwoPStatic) for a fixed
 *        pore network, sweeping the global capillary pressure through a full hysteresis cycle:
 *        up from `Problem.InitialPc` to `Problem.FinalPc` (drainage), then back down to
 *        `Problem.InitialPc` (imbibition). At every step of both legs, the resulting per-throat
 *        wetting-/non-wetting-phase transmissibilities are cached so they can later be imported
 *        into the single-phase upscaling solves (see main.cc), yielding relative-permeability
 *        curves that capture drainage/imbibition hysteresis.
 */
template<class TypeTag>
class PNMProblemStatic
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using StaticModel = Dumux::PoreNetwork::TwoPStatic<GridGeometry, Scalar>;

    using Transmissibility = Dumux::PoreNetwork::Detail::Transmissibility<PoreNetwork::TransmissibilityPatzekSilin<Scalar>,
                                                                          PoreNetwork::WettingLayerTransmissibility::RansohoffRadke<Scalar>,
                                                                          PoreNetwork::NonWettingPhaseTransmissibility::BakkeOren<Scalar>>;

public:

    //! The per-step results handed on to the single-phase upscaling solves.
    struct FlowPropertiesStatic
    {
        Scalar pcGlobal;
        Scalar averageSaturation;
        std::vector<std::array<Scalar, 2>> throatTransmissibility;
    };

    //! Total number of pc steps in the full drainage-then-imbibition cycle.
    int totalNumberOfSteps() const
    { return 2*numSteps_ + 1; }

    explicit PNMProblemStatic(const GridGeometry& gridGeometry, std::vector<bool> invasionState = {})
    : gridGeometry_(gridGeometry)
    , elementIsInvaded_(std::move(invasionState))
    {
        name_ = getParam<std::string>("Problem.Name");
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension");
        contactAngle_ = getParam<Scalar>("Problem.ContactAngle");
        inletPoreLabel_ = getParam<int>("Problem.InletPoreLabel");
        outletPoreLabel_ = getParam<int>("Problem.OutletPoreLabel");
        inletThroatLabel_ = inletPoreLabel_;
        outletThroatLabel_ = outletPoreLabel_;

        numSteps_ = getParam<int>("Problem.NumSteps");
        initialPc_ = getParam<Scalar>("Problem.InitialPc");
        finalPc_ = getParam<Scalar>("Problem.FinalPc");
        allowDraingeOfOutlet_ = getParam<bool>("Problem.AllowDraingeOfOutlet", true);

        flowPropertiesStatic_.resize(totalNumberOfSteps());
        const auto& leafGridView = gridGeometry_.gridView();

        if (elementIsInvaded_.empty())
            elementIsInvaded_.resize(leafGridView.size(0), false);
        elementIsTrapped_.resize(leafGridView.size(0), false);
        pcEntry_.resize(leafGridView.size(0));
        pcSnapOff_.resize(leafGridView.size(0));
        throatLabel_.resize(leafGridView.size(0));
        poreVolume_.resize(leafGridView.size(1));
        pc_.resize(leafGridView.size(1), 0.0);
        sw_.resize(leafGridView.size(1), 1.0);
        poreLabel_.resize(leafGridView.size(1));
        throatTransmissibility_.resize(leafGridView.size(0));

        const auto initialNumThroatsInvaded = static_cast<std::size_t>(
            std::count(elementIsInvaded_.begin(), elementIsInvaded_.end(), true));

        // helper function to evaluate the entry capillary pressure
        auto getPcEntry = [&](const std::size_t eIdx)
        {
            const Scalar throatRadius = gridGeometry_.throatInscribedRadius(eIdx);
            const auto shapeFactor = gridGeometry_.throatShapeFactor(eIdx);
            return PoreNetwork::ThresholdCapillaryPressures::pcEntry(surfaceTension_,
                                                                     contactAngle_,
                                                                     throatRadius,
                                                                     shapeFactor);
        };

        // helper function to evaluate the snap-off capillary pressure
        auto getPcSnapOff = [&](const std::size_t eIdx)
        {
            const Scalar throatRadius = gridGeometry_.throatInscribedRadius(eIdx);
            const auto shape = gridGeometry_.throatCrossSectionShape(eIdx);
            return PoreNetwork::ThresholdCapillaryPressures::pcSnapoff(surfaceTension_,
                                                                       contactAngle_,
                                                                       throatRadius,
                                                                       shape);
        };

        // fill the network properties
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);
            pcEntry_[eIdx] = getPcEntry(eIdx);
            pcSnapOff_[eIdx] = getPcSnapOff(eIdx);
            throatLabel_[eIdx] = gridGeometry_.throatLabel(eIdx);

            for (int i = 0; i < 2; ++i)
            {
                const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                poreLabel_[dofIdx] = gridGeometry_.poreLabel(dofIdx);
                poreVolume_[dofIdx] = gridGeometry_.poreVolume(dofIdx);
            }
        }

        totalPoreVolume_ = 0.0;
        for (std::size_t i = 0; i < poreVolume_.size(); ++i)
            if (poreLabel_[i] != inletPoreLabel_ && (poreLabel_[i] != outletPoreLabel_ || allowDraingeOfOutlet_))
                totalPoreVolume_ += poreVolume_[i];

        // prepare logfile
        logfileName_ = name_ + "_pc-s-curve.txt";
        logfile_.open(logfileName_);

        deltaPc_ = (finalPc_ - initialPc_) / numSteps_;
        pcGlobal_ = initialPc_;

        staticModel_ = std::make_unique<StaticModel>(gridGeometry_, pcEntry_, pcSnapOff_, throatLabel_,
                                                      inletThroatLabel_, outletThroatLabel_,
                                                      allowDraingeOfOutlet_, initialNumThroatsInvaded);
    }

    /*!
     * \brief Advances the invasion-percolation sweep by one step and caches the resulting
     *        per-throat phase transmissibilities and average saturation.
     */
    void pcSwStatic(int step)
    {
        std::cout << "Step " << step << ": Applying global pc of " << pcGlobal_ << " --> ";
        runStatic_();
        flowPropertiesStatic_[step] = {pcGlobal_, averageSaturation_, throatTransmissibility_};

        // step < numSteps_: still on the drainage leg; otherwise on the imbibition leg
        if (step < numSteps_)
            pcGlobal_ += deltaPc_;
        else
            pcGlobal_ -= deltaPc_;
    }

    /*!
     * \brief Plots the pc-S curve, with the drainage and imbibition legs shown as separate,
     *        distinctly styled and labeled curves so the hysteresis loop is visually clear.
     */
    void plot()
    {
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        gnuplot.addFileToPlot(logfileName_, "every ::0::" + std::to_string(numSteps_) +
                              " u 1:2 title 'Drainage' with linespoints lt 1 lc rgb '#1f77b4' pt 7 ps 0.6");
        gnuplot.addFileToPlot(logfileName_, "every ::" + std::to_string(numSteps_) +
                              " u 1:2 title 'Imbibition' with linespoints lt 2 lc rgb '#d62728' pt 5 ps 0.6");
        gnuplot.setXlabel("S_w [-]");
        gnuplot.setYlabel("p_c [Pa]");
        gnuplot.setXRange(0, 1.01);
        gnuplot.plot("plot");
    }

    auto sequenceWriter()
    {
        const auto& leafGridView = gridGeometry_.gridView();

        auto writer = std::make_shared<Dune::VTKWriter<GridView>>(leafGridView);
        Dune::VTKSequenceWriter<GridView> sequenceWriter(writer, "static");
        sequenceWriter.addCellData(pcEntry_, "pcEntry");
        sequenceWriter.addCellData(pcSnapOff_, "pcSnapOff");
        sequenceWriter.addCellData(elementIsInvaded_, "invaded");
        sequenceWriter.addCellData(elementIsTrapped_, "trapped");
        sequenceWriter.addCellData(throatLabel_, "throatLabel");
        sequenceWriter.addVertexData(poreLabel_, "poreLabel");
        sequenceWriter.addVertexData(pc_, "pc");
        sequenceWriter.addVertexData(sw_, "sw");
        sequenceWriter.addVertexData(poreVolume_, "poreVolume");

        return sequenceWriter;
    }

    const std::vector<std::array<Scalar, 2>>& throatTransmissibility() const
    { return throatTransmissibility_; }

    int numberOfSteps() const
    { return numSteps_; }

    const std::vector<FlowPropertiesStatic>& staticProperties() const
    { return flowPropertiesStatic_; }

    Scalar totalPoreVolume() const
    { return totalPoreVolume_; }

    const std::vector<bool>& invasionState() const
    { return elementIsInvaded_; }

private:

    //! Advances the invasion state by one step and updates the per-pore pc/Sw and the
    //! per-throat transmissibilities for the current global capillary pressure.
    void runStatic_()
    {
        const auto& leafGridView = gridGeometry_.gridView();

        // (a) update the throat invasion state for the current global capillary pressure
        staticModel_->updateInvasionState(elementIsInvaded_, elementIsTrapped_, pcGlobal_);

        // (b) update the per-pore capillary pressure/saturation, gated on the trapped state as
        // of the *previous* step (see dumux/porenetwork/2p/static/staticpnm.hh for why)
        std::fill(sw_.begin(), sw_.end(), 1.0);
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);
            if (elementIsInvaded_[eIdx] && !elementIsTrapped_[eIdx])
            {
                for (int i = 0; i < 2; ++i)
                {
                    const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                    pc_[dofIdx] = pcGlobal_;
                }
            }
        }

        averageSaturation_ = 0.0;
        auto fvGeometry = localView(gridGeometry_);
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bind(element);

            std::array<Scalar, 2> pcElement{};
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdx = scv.dofIndex();
                pcElement[scv.localDofIndex()] = pc_[dofIdx];

                if (pc_[dofIdx] > 0.0)
                {
                    using MaterialLaw = PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBodyDefault<PoreNetwork::Pore::Shape::cube>;
                    const Scalar poreRadius = gridGeometry_.poreInscribedRadius(dofIdx);

                    const auto params = MaterialLaw::BasicParams().setPoreInscribedRadius(poreRadius).setPoreShape(PoreNetwork::Pore::Shape::cube).setSurfaceTension(surfaceTension_);
                    auto fluidMatrixInteraction = makeFluidMatrixInteraction(MaterialLaw(params, MaterialLaw::RegularizationParams(), "SpatialParams"));
                    sw_[dofIdx] = fluidMatrixInteraction.sw(pc_[dofIdx]);
                }
                else
                    sw_[dofIdx] = 1.0;

                if (poreLabel_[dofIdx] != inletPoreLabel_ && (poreLabel_[dofIdx] != outletPoreLabel_ || allowDraingeOfOutlet_))
                    averageSaturation_ += scv.volume() * sw_[dofIdx];
            }

            calculateTransmissibility_(element, fvGeometry, pcElement);
        }
        averageSaturation_ /= totalPoreVolume_;

        // (c) recompute the trapped state for the next step, now that this step's saturation
        // update has already consumed the previous trapped state
        staticModel_->updateTrappedState(elementIsTrapped_, elementIsInvaded_);

        logfile_ << averageSaturation_ << " " << pcGlobal_ << std::endl;
        std::cout << staticModel_->numThroatsInvaded() << " of " << elementIsInvaded_.size()
                   << " throats invaded; Sw_avg " << averageSaturation_ << std::endl;
    }

    template<class Element, class FVElementGeometry>
    void calculateTransmissibility_(const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const std::array<Scalar, 2>& pc)
    {
        fluxVariablesCache_.update(gridGeometry_, element, pc);
        const auto eIdx = gridGeometry_.gridView().indexSet().index(element);
        const auto wPhaseIdx = fluxVariablesCache_.wPhaseIdx();
        const auto nPhaseIdx = fluxVariablesCache_.nPhaseIdx();

        const StaticUpscalingMockProblem problem{};
        const StaticUpscalingMockElementVolumeVariables elemVolVars{};
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (elementIsInvaded_[eIdx])
            {
                throatTransmissibility_[eIdx][wPhaseIdx] = Transmissibility::wettingLayerTransmissibility(element, fvGeometry, scvf, fluxVariablesCache_);
                throatTransmissibility_[eIdx][nPhaseIdx] = Transmissibility::nonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVariablesCache_);
            }
            else
            {
                throatTransmissibility_[eIdx][wPhaseIdx] = Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVariablesCache_, wPhaseIdx);
                // Not-yet-invaded throats carry no non-wetting phase, so their non-wetting
                // transmissibility should be (physically) zero. A literal zero would leave the
                // non-wetting-phase Darcy system singular wherever it disconnects part of the
                // network, so a small but nonzero placeholder is used instead. Keeping the
                // contrast to the real (invaded-throat) transmissibilities within about 6 orders
                // of magnitude -- rather than the 1e-10 used previously -- is enough to make the
                // non-wetting phase's contribution negligible while staying within what a direct
                // linear solver can resolve in double precision; 1e-10 was observed to push the
                // single-phase upscaling Newton solve to its precision floor and fail to converge
                // at low-saturation steps of the imbibition leg, where almost the whole network
                // is not yet (or no longer) invaded.
                static constexpr Scalar uninvadedNonWettingTransmissibilityFactor = 1e-6;
                throatTransmissibility_[eIdx][nPhaseIdx] = uninvadedNonWettingTransmissibilityFactor * throatTransmissibility_[eIdx][wPhaseIdx];
            }
        }
    }

    const GridGeometry& gridGeometry_;
    std::unique_ptr<StaticModel> staticModel_;
    FluxVariablesCache fluxVariablesCache_;
    std::ofstream logfile_;

    std::string name_;
    Scalar surfaceTension_;
    Scalar contactAngle_;
    int inletPoreLabel_;
    int outletPoreLabel_;
    int inletThroatLabel_;
    int outletThroatLabel_;

    int numSteps_;
    Scalar initialPc_;
    Scalar finalPc_;
    bool allowDraingeOfOutlet_;
    Scalar deltaPc_;
    Scalar pcGlobal_;

    std::vector<bool> elementIsInvaded_;
    std::vector<bool> elementIsTrapped_;
    std::vector<Scalar> pcEntry_;
    std::vector<Scalar> pcSnapOff_;
    std::vector<int> throatLabel_;
    std::vector<Scalar> poreVolume_;
    std::vector<Scalar> pc_;
    std::vector<Scalar> sw_;
    Scalar averageSaturation_;
    std::vector<int> poreLabel_;
    std::vector<std::array<Scalar, 2>> throatTransmissibility_;
    Scalar totalPoreVolume_;
    std::string logfileName_;
    std::vector<FlowPropertiesStatic> flowPropertiesStatic_;
};

} // end namespace Dumux::PoreNetwork

#endif
