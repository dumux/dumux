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

namespace Dumux::PoreNetwork::Detail {

template<class... TransmissibilityLawTypes>
struct Transmissibility : public TransmissibilityLawTypes... {};

}

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
        allowDraingeOfOutlet_ = getParam<bool>("Problem.AllowDraingeOfOutlet", false);

        pcRange_ = finalPc_ - initialPc_;
        deltaPc_ = pcRange_ / numSteps_;

        pcGlobal_ = initialPc_;


        const auto& leafGridView = gridGeometry_.gridView();
        elementIsInvaded_.resize(leafGridView.size(0), false);
        pcEntry_.resize(leafGridView.size(0));
        throatLabel_.resize(leafGridView.size(0));
        poreVolume_.resize(leafGridView.size(1));
        pc_.resize(leafGridView.size(1), 0.0);
        sw_.resize(leafGridView.size(1), 0.0);
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

    void pcSwStatic(DrainageModel& drainageModel)
    {
        const auto& leafGridView = gridGeometry_.gridView();

        fillNetworkProperties_(leafGridView);

        // // prepare logfile
        // std::ofstream logfile;
        // const auto logfileName = name_ + "_pc-s-curve.txt";
        // logfile.open(logfileName);

    //////////////////////////////////////////////
        // do the actual drainage process
            std::cout << ": Applying global pc of " << pcGlobal_ << " --> ";
            runDrainage_(leafGridView, drainageModel);
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
        static const auto name = getParam<std::string>("Problem.Name");
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
        Scalar averageSaturation = 0;
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
            const auto eIdx = leafGridView.indexSet().index(element);

            std::array<Scalar, 2> pcElement{};

            // get the saturation of each pore by inverting the local pore body pc-s curve
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdx = scv.dofIndex();

                pcElement[scv.localDofIndex()] = pc_[dofIdx];

                if (poreLabel_[dofIdx] == inletPoreLabel_ || (poreLabel_[dofIdx] == outletPoreLabel_ && !allowDraingeOfOutlet_))
                    continue;

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

                const Scalar partialPoreVolume = scv.volume();
                averageSaturation += partialPoreVolume*sw_[dofIdx];
            }
            fluxVariablesCache_.update(gridGeometry_, element, pcElement);
                // helper function to evaluate the transmissibility
            auto calculateTransmissibility = [&](std::size_t phaseIdx)
            {
                const auto problem = Dumux::Porenetwork::MockProblem{};
                const auto elemVolVars = Dumux::Porenetwork::MockElemVolVars{};
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    if (phaseIdx == fluxVariablesCache_.wPhaseIdx()) // wetting-wetting phase
                    {
                        return elementIsInvaded_[eIdx] ? Transmissibility::wettingLayerTransmissibility(element, fvGeometry, scvf, fluxVariablesCache_)
                                                        : Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVariablesCache_, phaseIdx);
                    }
                    else // non-wetting phase
                    {
                        return elementIsInvaded_[eIdx] ? Transmissibility::nonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVariablesCache_)
                                                        : 0.0;
                    }
                }
                return 0.0;
            };

            throatTransmissibility_[eIdx][0] = calculateTransmissibility(0);
            throatTransmissibility_[eIdx][1] = calculateTransmissibility(1);

            std::cout<<throatTransmissibility_[eIdx][0]<<"  "<<throatTransmissibility_[eIdx][1]<<std::endl;

        }
        averageSaturation /= totalPoreVolume_;

        // write to logfile
        logfile_ << averageSaturation << " " << pcGlobal_ << std::endl;

        // increment the global capillary pressure
        pcGlobal_ += deltaPc_;

        std::cout << drainageModel.numThroatsInvaded() << " of " << elementIsInvaded_.size() << " throats invaded; S_avg " << averageSaturation << std::endl;

    }

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
    std::vector<int> poreLabel_;
    std::vector<std::array<Scalar, 2>> throatTransmissibility_;
    Scalar totalPoreVolume_;

    std::string logfileName_;

};

}