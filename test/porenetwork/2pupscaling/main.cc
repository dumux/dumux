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

#include "properties.hh"

namespace Dumux{

template<class Scalar>
struct NetworkProperties{
    const Scalar surfaceTension = getParam<Scalar>("Problem.SurfaceTension");
    const Scalar contactAngle = getParam<Scalar>("Problem.ContactAngle");
    const int inletPoreLabel = getParam<int>("Problem.InletPoreLabel");
    const int outletPoreLabel = getParam<int>("Problem.OutletPoreLabel");
    const int inletThroatLabel = inletPoreLabel;
    const int outletThroatLabel = outletPoreLabel;
};

template<class Scalar>
struct DrainageData{
    const int numSteps = getParam<int>("Problem.NumSteps");
    const Scalar initialPc = getParam<Scalar>("Problem.InitialPc");
    const Scalar finalPc = getParam<Scalar>("Problem.FinalPc");
    const bool allowDraingeOfOutlet = getParam<bool>("Problem.AllowDraingeOfOutlet", false);

    const Scalar pcRange = finalPc - initialPc;
    const Scalar deltaPc = pcRange / numSteps;
};

template<class Scalar, class LeafGridView, class GridGeometry, class NetworkProperties>
struct SimulationData{

    SimulationData(const LeafGridView& leafGridView, const std::shared_ptr<GridGeometry> gridGeometry, const NetworkProperties& nP)
    :leafGridView_( leafGridView )
    {
        elementIsInvaded.resize(leafGridView_.size(0), false);
        pcEntry.resize(leafGridView_.size(0));
        throatLabel.resize(leafGridView_.size(0));
        poreVolume.resize(leafGridView_.size(1));
        pc.resize(leafGridView_.size(1), 0.0);
        sw.resize(leafGridView_.size(1), 0.0);
        poreLabel.resize(leafGridView_.size(1));
        throatTransmissibility.resize(leafGridView_.size(0));

        // helper function to evaluate the entry capillary pressure
        auto getPcEntry = [&](const std::size_t eIdx)
        {
            const Scalar throatRadius =  gridGeometry->throatInscribedRadius(eIdx);
            const auto shapeFactor = gridGeometry->throatShapeFactor(eIdx);
            return PoreNetwork::ThresholdCapillaryPressures::pcEntry(nP.surfaceTension,
                                                                     nP.contactAngle,
                                                                     throatRadius,
                                                                     shapeFactor);
        };

        // fill the network properties
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);
            pcEntry[eIdx] = getPcEntry(eIdx);
            throatLabel[eIdx] = gridGeometry->throatLabel(eIdx);

            for (int i = 0; i < 2; ++i)
            {
                const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                poreLabel[dofIdx] =  gridGeometry->poreLabel(dofIdx);
                poreVolume[dofIdx] =  gridGeometry->poreVolume(dofIdx);
            }
        }
    }
private:
    const LeafGridView& leafGridView_;

public:
    std::vector<bool> elementIsInvaded;
    std::vector<Scalar> pcEntry;
    std::vector<int> throatLabel;
    std::vector<Scalar> poreVolume;
    std::vector<Scalar> pc;
    std::vector<Scalar> sw;
    std::vector<int> poreLabel;
    std::vector<std::array<Scalar, 2>> throatTransmissibility;

};

template<class TypeTag>
void pcSw()
{
    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridManager = PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    // get some network properties
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    NetworkProperties<Scalar> nP;

    // define the drainge process
    DrainageData<Scalar> dData;

    // helper function to evaluate the entry capillary pressure
    auto getPcEntry = [&](const std::size_t eIdx)
    {
        const Scalar throatRadius =  gridGeometry->throatInscribedRadius(eIdx);
        const auto shapeFactor = gridGeometry->throatShapeFactor(eIdx);
        return PoreNetwork::ThresholdCapillaryPressures::pcEntry(nP.surfaceTension,
                                                                 nP.contactAngle,
                                                                 throatRadius,
                                                                 shapeFactor);
    };

    using Transmissibility = Dumux::PoreNetwork::Detail::Transmissibility<PoreNetwork::TransmissibilityPatzekSilin<Scalar>,
                                                                          PoreNetwork::WettingLayerTransmissibility::RansohoffRadke<Scalar>,
                                                                          PoreNetwork::NonWettingPhaseTransmissibility::BakkeOren<Scalar>>;


    // simulation data
    SimulationData<Scalar, decltype(leafGridView), GridGeometry, NetworkProperties<Scalar>> sData(leafGridView, gridGeometry, nP);

    // add vtk output
    static const auto name = getParam<std::string>("Problem.Name");
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    auto writer = std::make_shared<Dune::VTKWriter<GridView>>(leafGridView);
    Dune::VTKSequenceWriter<GridView> sequenceWriter(writer, name);
    sequenceWriter.addCellData(sData.pcEntry , "pcEntry");
    sequenceWriter.addCellData(sData.elementIsInvaded , "invaded");
    sequenceWriter.addCellData(sData.throatLabel , "throatLabel");
    sequenceWriter.addVertexData(sData.poreLabel , "poreLabel");
    sequenceWriter.addVertexData(sData.pc , "pc");
    sequenceWriter.addVertexData(sData.sw , "sw");
    sequenceWriter.addVertexData(sData.poreVolume , "poreVolume");

    // fill the network properties
    for (const auto& element : elements(leafGridView))
    {
        const auto eIdx = leafGridView.indexSet().index(element);
        sData.pcEntry[eIdx] = getPcEntry(eIdx);
        sData.throatLabel[eIdx] = gridGeometry->throatLabel(eIdx);

        for (int i = 0; i < 2; ++i)
        {
            const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
            sData.poreLabel[dofIdx] =  gridGeometry->poreLabel(dofIdx);
            sData.poreVolume[dofIdx] =  gridGeometry->poreVolume(dofIdx);
        }
    }

    // get the discrete capillary pressure steps
    Scalar pcGlobal = dData.initialPc;

    // get the drainage model
    PoreNetwork::TwoPStaticDrainage<GridGeometry, Scalar> drainageModel(*gridGeometry, sData.pcEntry, sData.throatLabel,
                                                                        nP.inletThroatLabel, nP.outletThroatLabel,
                                                                        dData.allowDraingeOfOutlet);

    // prepare logfile
    std::ofstream logfile;
    const auto logfileName = name + "_pc-s-curve.txt";
    logfile.open(logfileName);
    Scalar averageSaturation = 0;
    const Scalar totalPoreVolume = [&]()
    {
        Scalar result = 0.0;
        for (int i = 0; i < sData.poreVolume.size(); ++i)
            if (sData.poreLabel[i] != nP.inletPoreLabel && (sData.poreLabel[i] != nP.outletPoreLabel || dData.allowDraingeOfOutlet))
                result += sData.poreVolume[i];
        return result;
    }();

    std::cout << "total pore volume is " << totalPoreVolume << std::endl;

    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    FluxVariablesCache fluxVariablesCache;

    // do the actual drainage process
    for (int step = 0; step < dData.numSteps + 1; ++step)
    {
        std::cout << "Step " << step << ": Applying global pc of " << pcGlobal << " --> ";
        drainageModel.updateInvasionState(sData.elementIsInvaded, pcGlobal);

        // calculate the average saturation of the network
        averageSaturation = 0;
        auto fvGeometry = localView(*gridGeometry);
        std::fill(sData.sw.begin(), sData.sw.end(), 1.0);

        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);

            // update the capillary pressure of each pore body
            // pores connected to the invading phase receive the global capillary pressure, the otherones are set to zero
            if (sData.elementIsInvaded[eIdx])
            {
                for (int i = 0; i < 2; ++i)
                {
                    const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                    sData.pc[dofIdx] = pcGlobal;
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

                pcElement[scv.localDofIndex()] = sData.pc[dofIdx];

                if (sData.poreLabel[dofIdx] == nP.inletPoreLabel || (sData.poreLabel[dofIdx] == nP.outletPoreLabel && !dData.allowDraingeOfOutlet))
                    continue;

                if (sData.pc[dofIdx] > 0.0)
                {
                    using MaterialLaw = PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBodyDefault<PoreNetwork::Pore::Shape::cube>;
                    const Scalar poreRadius = gridGeometry->poreInscribedRadius(dofIdx);

                    const auto params = MaterialLaw::BasicParams().setPoreInscribedRadius(poreRadius).setPoreShape(PoreNetwork::Pore::Shape::cube).setSurfaceTension(nP.surfaceTension);
                    auto fluidMatrixInteraction = makeFluidMatrixInteraction(MaterialLaw(params, MaterialLaw::RegularizationParams(), "SpatialParams"));
                    sData.sw[dofIdx] = fluidMatrixInteraction.sw(sData.pc[dofIdx]);
                }
                else
                    sData.sw[dofIdx] = 1.0;

                const Scalar partialPoreVolume = scv.volume();
                averageSaturation += partialPoreVolume*sData.sw[dofIdx];
            }
            fluxVariablesCache.update(gridGeometry, element, pcElement);
                // helper function to evaluate the transmissibility
            auto calculateTransmissibility = [&](std::size_t phaseIdx)
            {
                const auto problem = Dumux::Porenetwork::MockProblem{};
                const auto elemVolVars = Dumux::Porenetwork::MockElemVolVars{};
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    if (phaseIdx == fluxVariablesCache.wPhaseIdx()) // wetting-wetting phase
                    {
                        return sData.elementIsInvaded[eIdx] ? Transmissibility::wettingLayerTransmissibility(element, fvGeometry, scvf, fluxVariablesCache)
                                                            : Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVariablesCache, phaseIdx);
                    }
                    else // non-wetting phase
                    {
                        return sData.elementIsInvaded[eIdx] ? Transmissibility::nonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVariablesCache)
                                                            : 0.0;
                    }
                }
                return 0.0;
            };

            sData.throatTransmissibility[eIdx][0] = calculateTransmissibility(0);
            sData.throatTransmissibility[eIdx][1] = calculateTransmissibility(1);

            std::cout<<sData.throatTransmissibility[eIdx][0]<<"  "<<sData.throatTransmissibility[eIdx][1]<<std::endl;

        }
        averageSaturation /= totalPoreVolume;

        // write to logfile
        logfile << averageSaturation << " " << pcGlobal << std::endl;

        // increment the global capillary pressure
        pcGlobal += dData.deltaPc;

        std::cout << drainageModel.numThroatsInvaded() << " of " << sData.elementIsInvaded.size() << " throats invaded; S_avg " << averageSaturation << std::endl;
        sequenceWriter.write(step);
    }

    //plot the pc-S curve, if desired
#ifdef HAVE_GNUPLOT
    if (getParam<bool>("Problem.PlotPcS"))
    {
        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        gnuplot.addFileToPlot(logfileName);
        gnuplot.setXlabel("S_w [-]");
        gnuplot.setYlabel("p_c [P]");
        gnuplot.setXRange(0, 1.01);
        gnuplot.plot("plot");
    }
#endif

};

}

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::DrainageProblem;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    pcSw<TypeTag>();

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    // timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }
}
