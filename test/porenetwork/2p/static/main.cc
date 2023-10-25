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

#include <dumux/porenetwork/2p/static/staticdrainge.hh>
#include <dumux/io/gnuplotinterface.hh>

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DrainageProblem { using InheritsFrom = std::tuple<GridProperties, ModelProperties>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::DrainageProblem>
{
private:
    static constexpr bool enableCache = false;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::PoreNetwork::GridGeometry<Scalar, GridView, enableCache>;
};

} // end namespace Dumux::Properties


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
    const Scalar surfaceTension = getParam<Scalar>("Problem.SurfaceTension");
    const Scalar contactAngle = getParam<Scalar>("Problem.ContactAngle");
    const int inletPoreLabel = getParam<int>("Problem.InletPoreLabel");
    const int outletPoreLabel = getParam<int>("Problem.OutletPoreLabel");
    const int inletThroatLabel = inletPoreLabel;
    const int outletThroatLabel = outletPoreLabel;

    // define the drainge process
    const int numSteps = getParam<int>("Problem.NumSteps");
    const Scalar initialPc = getParam<Scalar>("Problem.InitialPc");
    const Scalar finalPc = getParam<Scalar>("Problem.FinalPc");
    const bool allowDraingeOfOutlet = getParam<bool>("Problem.AllowDraingeOfOutlet", false);

    // helper function to evaluate the entry capillary pressure
    auto getPcEntry = [&](const std::size_t eIdx)
    {
        const Scalar throatRadius =  gridGeometry->throatInscribedRadius(eIdx);
        const auto shapeFactor = gridGeometry->throatShapeFactor(eIdx);
        return PoreNetwork::ThresholdCapillaryPressures::pcEntry(surfaceTension,
                                                                 contactAngle,
                                                                 throatRadius,
                                                                 shapeFactor);
    };

    // simulation data
    std::vector<bool> elementIsInvaded(leafGridView.size(0), false);
    std::vector<Scalar> pcEntry(leafGridView.size(0));
    std::vector<int> throatLabel(leafGridView.size(0));
    std::vector<Scalar> poreVolume(leafGridView.size(1));
    std::vector<Scalar> pc(leafGridView.size(1), 0.0);
    std::vector<Scalar> sw(leafGridView.size(1), 0.0);
    std::vector<int> poreLabel(leafGridView.size(1));

    // add vtk output
    static const auto name = getParam<std::string>("Problem.Name");
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    auto writer = std::make_shared<Dune::VTKWriter<GridView>>(leafGridView);
    Dune::VTKSequenceWriter<GridView> sequenceWriter(writer, name);
    sequenceWriter.addCellData(pcEntry , "pcEntry");
    sequenceWriter.addCellData(elementIsInvaded , "invaded");
    sequenceWriter.addCellData(throatLabel , "throatLabel");
    sequenceWriter.addVertexData(poreLabel , "poreLabel");
    sequenceWriter.addVertexData(pc , "pc");
    sequenceWriter.addVertexData(sw , "sw");
    sequenceWriter.addVertexData(poreVolume , "poreVolume");

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

    // get the discrete capillary pressure steps
    const Scalar pcRange = finalPc - initialPc;
    const Scalar deltaPc = pcRange / numSteps;
    Scalar pcGlobal = initialPc;

    // get the drainage model
    PoreNetwork::TwoPStaticDrainage<GridGeometry, Scalar> drainageModel(*gridGeometry, pcEntry, throatLabel,
                                                                        inletThroatLabel, outletThroatLabel,
                                                                        allowDraingeOfOutlet);

    // prepare logfile
    std::ofstream logfile;
    const auto logfileName = name + "_pc-s-curve.txt";
    logfile.open(logfileName);
    Scalar averageSaturation = 0;
    const Scalar totalPoreVolume = [&]()
    {
        Scalar result = 0.0;
        for (int i = 0; i < poreVolume.size(); ++i)
            if (poreLabel[i] != inletPoreLabel && (poreLabel[i] != outletPoreLabel || allowDraingeOfOutlet))
                result += poreVolume[i];
        return result;
    }();

    std::cout << "total pore volume is " << totalPoreVolume << std::endl;

    // do the actual drainage process
    for (int step = 0; step < numSteps + 1; ++step)
    {
        std::cout << "Step " << step << ": Applying global pc of " << pcGlobal << " --> ";
        drainageModel.updateInvasionState(elementIsInvaded, pcGlobal);

        // calculate the average saturation of the network
        averageSaturation = 0;
        auto fvGeometry = localView(*gridGeometry);
        std::fill(sw.begin(), sw.end(), 1.0);

        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);

            // update the capillary pressure of each pore body
            // pores connected to the invading phase receive the global capillary pressure, the otherones are set to zero
            if (elementIsInvaded[eIdx])
            {
                for (int i = 0; i < 2; ++i)
                {
                    const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
                    pc[dofIdx] = pcGlobal;
                }
            }
        }

        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bind(element);

            // get the saturation of each pore by inverting the local pore body pc-s curve
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdx = scv.dofIndex();
                if (poreLabel[dofIdx] == inletPoreLabel || (poreLabel[dofIdx] == outletPoreLabel && !allowDraingeOfOutlet))
                    continue;

                if (pc[dofIdx] > 0.0)
                {
                    using MaterialLaw = PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBodyDefault<PoreNetwork::Pore::Shape::cube>;
                    const Scalar poreRadius = gridGeometry->poreInscribedRadius(dofIdx);

                    const auto params = MaterialLaw::BasicParams().setPoreInscribedRadius(poreRadius).setPoreShape(PoreNetwork::Pore::Shape::cube).setSurfaceTension(surfaceTension);
                    auto fluidMatrixInteraction = makeFluidMatrixInteraction(MaterialLaw(params, MaterialLaw::RegularizationParams(), "SpatialParams"));
                    sw[dofIdx] = fluidMatrixInteraction.sw(pc[dofIdx]);
                }
                else
                    sw[dofIdx] = 1.0;

                const Scalar partialPoreVolume = scv.volume();
                averageSaturation += partialPoreVolume*sw[dofIdx];
            }
        }
        averageSaturation /= totalPoreVolume;

        // write to logfile
        logfile << averageSaturation << " " << pcGlobal << std::endl;

        // increment the global capillary pressure
        pcGlobal += deltaPc;

        std::cout << drainageModel.numThroatsInvaded() << " of " << elementIsInvaded.size() << " throats invaded; S_avg " << averageSaturation << std::endl;
        sequenceWriter.write(step);
    }

    //plot the pc-S curve, if desired
#ifdef DUMUX_HAVE_GNUPLOT
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
