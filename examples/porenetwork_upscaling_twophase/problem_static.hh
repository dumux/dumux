// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Problem file for upscaling using the pore network model
*/
#include <config.h>
#include <ctime>
#include <iostream>
#include<algorithm>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/common/properties/model.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/discretization/porenetwork/gridgeometry.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

#include <dumux/porenetwork/2p/static/staticpnm.hh>
#include <dumux/io/gnuplotinterface.hh>

namespace Dumux::Porenetwork{

    struct MockProblem{};
    struct MockElemVolVars{};

    struct MockFVElementGeometry
    {
        struct SubControlVolumeFace{};
    };
}

namespace Dumux::PoreNetwork{

template<class TypeTag>
class PNMProblemStatic{

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using StaticModel = Dumux::PoreNetwork::TwoPStatic<GridGeometry, Scalar, DRAINAGE>;


    using Transmissibility = Dumux::PoreNetwork::Detail::Transmissibility<PoreNetwork::TransmissibilityPatzekSilin<Scalar>,
                                                                          PoreNetwork::WettingLayerTransmissibility::RansohoffRadke<Scalar>,
                                                                          PoreNetwork::NonWettingPhaseTransmissibility::BakkeOren<Scalar>>;

    struct FlowPropertiesStatic
    {
        Scalar pcGlobal;
        Scalar averageSaturation;
        std::vector<std::array<Scalar, 2>> throatTransmissibility;
    };

public:

    PNMProblemStatic(const GridGeometry& gridGeometry, bool isDrainage, std::vector<bool> invasionState = {})
    :gridGeometry_(gridGeometry)
    , elementIsInvaded_(invasionState)
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
        isDrainageProcess_ = isDrainage;
        allowDraingeOfOutlet_ = getParam<bool>("Problem.AllowDraingeOfOutlet", true);

        flowPropertiesStatic_.resize(numSteps_+1);
        const auto& leafGridView = gridGeometry_.gridView();

        if (elementIsInvaded_.empty())
            elementIsInvaded_.resize(leafGridView.size(0), !isDrainageProcess_);
        pcEntry_.resize(leafGridView.size(0));
        pcSnapOff_.resize(leafGridView.size(0));
        throatLabel_.resize(leafGridView.size(0));
        poreVolume_.resize(leafGridView.size(1));
        pc_.resize(leafGridView.size(1), 0.0);
        sw_.resize(leafGridView.size(1), 1.0);
        poreLabel_.resize(leafGridView.size(1));
        throatTransmissibility_.resize(leafGridView.size(0));
        throatVolume_.resize(leafGridView.size(0));

        initialNumThroatsInvaded_ = std::count_if(elementIsInvaded_.begin(), elementIsInvaded_.end(), [](bool value) { return value; });

        // helper function to evaluate the entry capillary pressure
        auto getPcEntry = [&](const std::size_t eIdx)
        {
            const Scalar throatRadius =  gridGeometry_.throatInscribedRadius(eIdx);
            const auto shapeFactor = gridGeometry_.throatShapeFactor(eIdx);
            return PoreNetwork::ThresholdCapillaryPressures::pcEntry(surfaceTension_,
                                                                     contactAngle_,
                                                                     throatRadius,
                                                                     shapeFactor);
        };

        // helper function to evaluate the entry capillary pressure
        auto getPcSnapOff = [&](const std::size_t eIdx)
        {
            const Scalar throatRadius =  gridGeometry_.throatInscribedRadius(eIdx);
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
            const Scalar throatRadius =  gridGeometry_.throatInscribedRadius(eIdx);
            const Scalar throatLength =  gridGeometry_.throatLength(eIdx);
            throatVolume_[eIdx] = throatRadius * throatRadius * throatLength; //square cross-section

            for (int i = 0; i < 2; ++i)
            {
                const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                poreLabel_[dofIdx] =  gridGeometry_.poreLabel(dofIdx);
                poreVolume_[dofIdx] =  gridGeometry_.poreVolume(dofIdx);
            }
        }
        totalPoreVolume_ = 0.0;
        for (int i = 0; i < poreVolume_.size(); ++i)
            if (poreLabel_[i] != inletPoreLabel_ && (poreLabel_[i] != outletPoreLabel_ || allowDraingeOfOutlet_))
                totalPoreVolume_ += poreVolume_[i];

        totalThroatVolume_ = std::accumulate(throatVolume_.begin(), throatVolume_.end(), 0.0);

        // prepare logfile
        logfileName_ = name_ + "_pc-s-curve.txt";
        logfile_.open(logfileName_);

        // set the final to zero and the initial to the maximum pc-entry for imbibition process
        if (!isDrainageProcess_)
        {
            initialPc_ = getParam<Scalar>("Problem.FinalPc");
            finalPc_ = getParam<Scalar>("Problem.InitialPc");
        }

        pcRange_ = finalPc_ - initialPc_;
        deltaPc_ = pcRange_ / numSteps_;
        pcGlobal_ = initialPc_;

    }


    template<class StaticModel>
    void pcSwStatic(StaticModel& staticModel, int step)
    {
        const auto& leafGridView = gridGeometry_.gridView();
        Scalar pcTemp = pcGlobal_;

        flowPropertiesStatic_[step].throatTransmissibility.resize(leafGridView.size(0));
        std::cout << "Step " << step << " --> "<<": Applying global pc of " << pcGlobal_ << " --> ";
        runStatic_(leafGridView, staticModel);
        flowPropertiesStatic_[step] = {pcGlobal_, averageSaturation_, throatTransmissibility_};

        // increment the global capillary pressure
        pcTemp += deltaPc_;

        pcGlobal_ = std::max(pcTemp, 0.0);
    }

    //plot the pc-S curve, if desired
    void plot()
    {
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        gnuplot.addFileToPlot(logfileName_);
        gnuplot.setXlabel("S_w [-]");
        gnuplot.setYlabel("p_c [P]");
        gnuplot.setXRange(0, 1.01);
        gnuplot.plot("plot");
    }

    auto sequenceWriter()
    {
        const auto& leafGridView = gridGeometry_.gridView();

        // add vtk output
        static const auto name = "static";
        auto writer = std::make_shared<Dune::VTKWriter<GridView>>(leafGridView);
        Dune::VTKSequenceWriter<GridView> sequenceWriter(writer, name);
        sequenceWriterAddCell_(sequenceWriter);

        return sequenceWriter;
    }

    auto drainageModelStatic()
    {

        // get the drainage model
        return DrainageModel(gridGeometry_, pcEntry_, throatLabel_,
                             inletThroatLabel_, outletThroatLabel_,
                             allowDraingeOfOutlet_);
    }

    auto imbibitionModelStatic()
    {

        // get the imbibition model
        return ImbibitionModel(gridGeometry_, pcEntry_, pcSnapOff_, throatLabel_,
                               inletThroatLabel_, outletThroatLabel_,
                               allowDraingeOfOutlet_);
    }

    template<bool IsDrainage>
    auto staticModel()
    {
        using StaticModel = Dumux::PoreNetwork::TwoPStatic<GridGeometry, Scalar, IsDrainage>;
        return StaticModel(gridGeometry_, pcEntry_, pcSnapOff_, throatLabel_,
                               inletThroatLabel_, outletThroatLabel_, initialNumThroatsInvaded_,
                               allowDraingeOfOutlet_);
    }

    const auto throatTransmissibility() const
    { return throatTransmissibility_; }

    auto numberOfSteps()
    { return numSteps_; };

    auto staticProperties()
    { return flowPropertiesStatic_; };

    Scalar totalPoreVolume()
    { return totalPoreVolume_; }

    bool isDrainage()
    { return isDrainageProcess_;}

    auto invasionState()
    { return elementIsInvaded_; }

private:

    void sequenceWriterAddCell_(Dune::VTKSequenceWriter<GridView>& sequenceWriter)
    {
        sequenceWriter.addCellData(pcEntry_ , "pcEntry");
        sequenceWriter.addCellData(pcSnapOff_ , "pcSnapOff");
        sequenceWriter.addCellData(elementIsInvaded_ , "invaded");
        sequenceWriter.addCellData(throatLabel_ , "throatLabel");
        sequenceWriter.addVertexData(poreLabel_ , "poreLabel");
        sequenceWriter.addVertexData(pc_ , "pc");
        sequenceWriter.addVertexData(sw_ , "sw");
        sequenceWriter.addVertexData(poreVolume_ , "poreVolume");
    }

    template<class LeafGridView>
    void fillNetworkProperties_(const LeafGridView& leafGridView)
    {
        // helper function to evaluate the entry capillary pressure
        auto getPcEntry = [&](const std::size_t eIdx)
        {
            const Scalar throatRadius =  gridGeometry_.throatInscribedRadius(eIdx);
            const auto shapeFactor = gridGeometry_.throatShapeFactor(eIdx);
            return PoreNetwork::ThresholdCapillaryPressures::pcEntry(surfaceTension_,
                                                                     contactAngle_,
                                                                     throatRadius,
                                                                     shapeFactor);
        };

        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);
            pcEntry_[eIdx] = getPcEntry(eIdx);
            throatLabel_[eIdx] = gridGeometry_.throatLabel(eIdx);

            for (int i = 0; i < 2; ++i)
            {
                const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                poreLabel_[dofIdx] =  gridGeometry_.poreLabel(dofIdx);
                poreVolume_[dofIdx] =  gridGeometry_.poreVolume(dofIdx);
            }
        }
    }

    template<class LeafGridView, class StaticModel>
    void runStatic_(const LeafGridView& leafGridView, StaticModel& staticModel)
    {
        std::fill(sw_.begin(), sw_.end(), 0.0);

        staticModel.updateInvasionState(elementIsInvaded_, pcGlobal_);

        // calculate the average saturation of the network
        averageSaturation_ = 0;
        auto fvGeometry = localView(gridGeometry_);

        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);

            // update the capillary pressure of each pore body
            // pores connected to the invading phase receive the global capillary pressure, the otherones are set to zero
            if (elementIsInvaded_[eIdx])
            {
                for (int i = 0; i < 2; ++i)
                {
                    const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                    pc_[dofIdx] = pcGlobal_;
                }
            }
        }

        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bind(element);

            std::array<Scalar, 2> pcElement{};

            // get the saturation of each pore by inverting the local pore body pc-s curve
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
                {
                    const Scalar partialPoreVolume = scv.volume();
                    averageSaturation_ += partialPoreVolume*sw_[dofIdx];
                }
            }

            calculateTransmissibility_(element, fvGeometry, pcElement);

        }
        averageSaturation_ /= totalPoreVolume_;

        // write to logfile
        logfile_ << averageSaturation_ << " " << pcGlobal_ << std::endl;

        std::cout << staticModel.numThroatsInvaded() << " of " << elementIsInvaded_.size() << " throats invaded; Sw_avg " << averageSaturation_ << std::endl;

    }

    template<class Element, class FVElementGeometry>
    void calculateTransmissibility_(const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const std::array<Scalar,2> pc)
    {
        fluxVariablesCache_.update(gridGeometry_, element, pc);
        const auto eIdx = gridGeometry_.gridView().indexSet().index(element);
        const auto wPhaseIdx = fluxVariablesCache_.wPhaseIdx();
        const auto nPhaseIdx = fluxVariablesCache_.nPhaseIdx();

        const auto problem = Dumux::Porenetwork::MockProblem{};
        const auto elemVolVars = Dumux::Porenetwork::MockElemVolVars{};
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
                throatTransmissibility_[eIdx][nPhaseIdx] = 1e-10 * throatTransmissibility_[eIdx][wPhaseIdx];
            }
        }
    };

    const GridGeometry& gridGeometry_;
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
    bool isDrainageProcess_;
    bool allowDraingeOfOutlet_;
    std::size_t initialNumThroatsInvaded_;
    Scalar pcRange_;
    Scalar deltaPc_;
    Scalar pcGlobal_;

    std::vector<bool> elementIsInvaded_;
    std::vector<Scalar> pcEntry_;
    std::vector<Scalar> pcSnapOff_;
    std::vector<int> throatLabel_;
    std::vector<Scalar> poreVolume_;
    std::vector<Scalar> throatVolume_;
    std::vector<Scalar> pc_;
    std::vector<Scalar> sw_;
    Scalar averageSaturation_;
    std::vector<int> poreLabel_;
    std::vector<std::array<Scalar, 2>> throatTransmissibility_;
    Scalar totalPoreVolume_;
    Scalar totalThroatVolume_;
    std::string logfileName_;
    std::vector<FlowPropertiesStatic> flowPropertiesStatic_;

};

}
