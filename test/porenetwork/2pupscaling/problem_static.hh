// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the pore network model
*/
#include <config.h>
#include <ctime>
#include <iostream>
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

#include <dumux/porenetwork/2p/static/staticdrainge.hh>
#include <dumux/io/gnuplotinterface.hh>

// namespace Dumux::PoreNetwork::Detail {

// template<class... TransmissibilityLawTypes>
// struct Transmissibility : public TransmissibilityLawTypes... {};

// }

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
class DrainageProblemStatic{

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;

    using DrainageModel = Dumux::PoreNetwork::TwoPStaticDrainage<GridGeometry, Scalar>;

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

    DrainageProblemStatic(const GridGeometry& gridGeometry)
    :gridGeometry_(gridGeometry)
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

        pcRange_ = finalPc_ - initialPc_;
        deltaPc_ = pcRange_ / numSteps_;

        pcGlobal_ = initialPc_;

        flowPropertiesStatic_.resize(numSteps_+1);
        const auto& leafGridView = gridGeometry_.gridView();
        elementIsInvaded_.resize(leafGridView.size(0), false);
        pcEntry_.resize(leafGridView.size(0));
        throatLabel_.resize(leafGridView.size(0));
        poreVolume_.resize(leafGridView.size(1));
        pc_.resize(leafGridView.size(1), 0.0);
        sw_.resize(leafGridView.size(1), 1.0);
        poreLabel_.resize(leafGridView.size(1));
        throatTransmissibility_.resize(leafGridView.size(0));

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

        // fill the network properties
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
        totalPoreVolume_ = 0.0;
        for (int i = 0; i < poreVolume_.size(); ++i)
            if (poreLabel_[i] != inletPoreLabel_ && (poreLabel_[i] != outletPoreLabel_ || allowDraingeOfOutlet_))
                totalPoreVolume_ += poreVolume_[i];


        // prepare logfile
        logfileName_ = name_ + "_pc-s-curve.txt";
        logfile_.open(logfileName_);

    }

    void pcSwStatic(DrainageModel& drainageModel, int step)
    {
        const auto& leafGridView = gridGeometry_.gridView();

        fillNetworkProperties_(leafGridView);

        // // prepare logfile
        // std::ofstream logfile;
        // const auto logfileName = name_ + "_pc-s-curve.txt";
        // logfile.open(logfileName);

        // for (int step = 0; step < numSteps_ + 1; ++step)
        // {
        flowPropertiesStatic_[step].throatTransmissibility.resize(leafGridView.size(0));
        std::cout << "Step " << step << " --> "<<": Applying global pc of " << pcGlobal_ << " --> ";
        runDrainage_(leafGridView, drainageModel);
        //sequenceWriter.write(step);
        flowPropertiesStatic_[step] = {pcGlobal_, averageSaturation_, throatTransmissibility_};

        // increment the global capillary pressure
        pcGlobal_ += deltaPc_;
        //}

        // return flowPropertiesStatic_;
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

        // add vtk output TODO: Bring it outside the loop
        static const auto name = "static";//getParam<std::string>("Problem.Name");
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

    const auto throatTransmissibility() const
    { return throatTransmissibility_; }

    auto numberOfSteps()
    { return numSteps_; };

    auto staticProperties()
    { return flowPropertiesStatic_; };

    Scalar totalPoreVolume()
    { return totalPoreVolume_; }
private:

    void sequenceWriterAddCell_(Dune::VTKSequenceWriter<GridView>& sequenceWriter)
    {
        sequenceWriter.addCellData(pcEntry_ , "pcEntry");
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

    template<class LeafGridView>
    void runDrainage_(const LeafGridView& leafGridView, DrainageModel& drainageModel)
    {
        drainageModel.updateInvasionState(elementIsInvaded_, pcGlobal_);

        // calculate the average saturation of the network
        averageSaturation_ = 0;
        auto fvGeometry = localView(gridGeometry_);
        std::fill(sw_.begin(), sw_.end(), 1.0);

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

                // if (poreLabel_[dofIdx] == inletPoreLabel_ || (poreLabel_[dofIdx] == outletPoreLabel_ && !allowDraingeOfOutlet_))
                //     continue;

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

        std::cout << drainageModel.numThroatsInvaded() << " of " << elementIsInvaded_.size() << " throats invaded; S_avg " << averageSaturation_ << std::endl;

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
    bool allowDraingeOfOutlet_;
    Scalar pcRange_;
    Scalar deltaPc_;
    Scalar pcGlobal_;

    std::vector<bool> elementIsInvaded_;
    std::vector<Scalar> pcEntry_;
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

}
