// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// ## The main program (`main.cc`)
// This file contains the main program flow. In this example, we use a quasi-static
// two-phase Pore-Network-Model to determine the capillary pressure-saturation (pc-Sw)
// curve of a pore network via invasion percolation with snap-off and non-wetting-phase
// trapping.
//
// Unlike most DuMu<sup>x</sup> examples, this one does not assemble and solve a PDE:
// `Dumux::PoreNetwork::TwoPStatic` (see
// [`staticpnm.hh`](../../dumux/porenetwork/2p/static/staticpnm.hh)) performs a purely
// topological invasion-percolation sweep for a given global capillary pressure, so the
// main function is a straight-line driver rather than the usual grid/problem/assembler/
// solver machinery. At a glance, `main()` does the following:
//
// 1. build the pore-network grid and grid geometry from `params.input`,
// 2. precompute each throat's entry and snap-off capillary pressures once, since they only depend on the fixed network geometry, not on the current capillary pressure,
// 3. step a global capillary pressure `pcGlobal` from `InitialPc` up to `FinalPc` (drainage) and back down to `InitialPc` (imbibition); at each step it (a) updates the throat invasion state for `pcGlobal` via `updateInvasionState`, (b) updates each pore's local capillary pressure/saturation, gated on the *previous* step's trapped state, and (c) recomputes the trapped state for the next step via `updateTrappedState`,
// 4. writes the pore-volume-averaged saturation and `pcGlobal` of every step to `<name>_pc-s-curve.txt` (and to a `.vtp`/`.pvd` sequence for ParaView).
//
// See the [equations and pseudocode in the main description](../README.md#mathematical-and-numerical-model)
// for the physics and algorithm behind steps 2 and 3.
// [[content]]
// ### Includes
// [[details]] includes
// [[codeblock]]
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/discretization/porenetwork/gridgeometry.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/gnuplotinterface.hh>

#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>

#include <dumux/porenetwork/2p/static/staticpnm.hh>
// [[/codeblock]]
// [[/details]]
//
// ### Compile-time settings
// We create a new type tag for our simulation and inherit the properties needed
// for a pore-network grid. `TwoPStatic` (used further down) is not a full-fledged
// DuMu<sup>x</sup> model -- it does not need `Problem`/`GridVariables`/assembler
// property specializations -- so only the `Grid` and `GridGeometry` properties
// are set here.
// [[codeblock]]
namespace Dumux::Properties {

namespace TTag {
struct PNMTwoPStatic { using InheritsFrom = std::tuple<GridProperties, ModelProperties>; };
} // end namespace TTag

// the throats (elements) of the pore-network grid live on a 1d network embedded in 3d space
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMTwoPStatic> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PNMTwoPStatic>
{
private:
    static constexpr bool enableCache = false;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::PoreNetwork::GridGeometry<Scalar, GridView, enableCache>;
};

} // end namespace Dumux::Properties
// [[/codeblock]]
//
// ### Inverting a pore's local pc-Sw law for its actual shape
// The pore-local pc-Sw closed-form relation (`TwoPLocalRulesPlatonicBodyDefault`) is only
// implemented for "platonic body" pore shapes (tetrahedron, cube, octahedron, dodecahedron,
// icosahedron), see `dumux/material/fluidmatrixinteractions/porenetwork/pore/2p`. Since that
// relation is a class template parameterized on the shape at compile time, while the shape
// actually assigned to a pore (`Grid.PoreGeometry` in the input file) is only known at
// runtime, we need a small runtime-to-compile-time dispatcher: it instantiates the matching
// template for each supported shape and throws `Dune::NotImplemented` for any shape without
// an implemented pc-Sw relation (e.g. Sphere, Circle, Cylinder), rather than silently
// defaulting to some other shape's curve.
// [[codeblock]]
namespace {

template<class Scalar>
Scalar poreSaturation(Dumux::PoreNetwork::Pore::Shape shape,
                      const Scalar poreRadius, const Scalar surfaceTension, const Scalar pc)
{
    using namespace Dumux;

    auto invert = [&](auto shapeTag)
    {
        constexpr auto s = decltype(shapeTag)::value;
        using MaterialLaw = PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBodyDefault<s>;
        using BasicParams = typename MaterialLaw::BasicParams;
        using RegularizationParams = typename MaterialLaw::RegularizationParams;

        const auto params = BasicParams().setPoreInscribedRadius(poreRadius).setPoreShape(s).setSurfaceTension(surfaceTension);
        auto fluidMatrixInteraction = makeFluidMatrixInteraction(MaterialLaw(params, RegularizationParams(), "SpatialParams"));
        return fluidMatrixInteraction.sw(pc);
    };

    using PoreNetwork::Pore::Shape;
    switch (shape)
    {
        case Shape::tetrahedron:  return invert(std::integral_constant<Shape, Shape::tetrahedron>{});
        case Shape::cube:         return invert(std::integral_constant<Shape, Shape::cube>{});
        case Shape::octahedron:   return invert(std::integral_constant<Shape, Shape::octahedron>{});
        case Shape::dodecahedron: return invert(std::integral_constant<Shape, Shape::dodecahedron>{});
        case Shape::icosahedron:  return invert(std::integral_constant<Shape, Shape::icosahedron>{});
        default:
            DUNE_THROW(Dune::NotImplemented,
                       "No pc-Sw relation available for pore shape '" << PoreNetwork::Pore::shapeToString(shape)
                       << "'. Supported shapes: Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron "
                       << "(see dumux/material/fluidmatrixinteractions/porenetwork/pore/2p).");
    }
}

} // end anonymous namespace
// [[/codeblock]]
//
// ### The main function
// [[details]] main
// [[codeblock]]
int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::PNMTwoPStatic;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and the input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // ### Create the grid and the grid geometry
    // The grid is a randomly generated pore network read from the parameters in the
    // `[Grid]` group of `params.input` (the same network-generation parameters as in the
    // [pore-network upscaling example](../porenetwork_upscaling/README.md), so both
    // examples operate on the same kind of network).
    // [[codeblock]]
    using GridManager = PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);
    // [[/codeblock]]

    // ### Read fluid and process parameters
    // [[codeblock]]
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const Scalar surfaceTension = getParam<Scalar>("Problem.SurfaceTension");
    const Scalar contactAngle = getParam<Scalar>("Problem.ContactAngle");
    const int inletPoreLabel = getParam<int>("Problem.InletPoreLabel");
    const int outletPoreLabel = getParam<int>("Problem.OutletPoreLabel");
    const int inletThroatLabel = inletPoreLabel;
    const int outletThroatLabel = outletPoreLabel;

    // the pc-S curve is traced by numSteps discrete pc increments from InitialPc
    // up to FinalPc (drainage), followed by the same number of decrements back
    // down to InitialPc (imbibition)
    const int numSteps = getParam<int>("Problem.NumSteps");
    const Scalar initialPc = getParam<Scalar>("Problem.InitialPc");
    const Scalar finalPc = getParam<Scalar>("Problem.FinalPc");
    const bool allowDraingeOfOutlet = getParam<bool>("Problem.AllowDraingeOfOutlet", false);
    // [[/codeblock]]

    // ### Precompute the throat entry and snap-off capillary pressures
    // These are the two threshold pressures that drive the invasion-percolation
    // algorithm: `pcEntry` is the Young-Laplace/MS-P threshold at which the
    // non-wetting phase can invade a throat (drainage), `pcSnapOff` is the Roof
    // snap-off threshold below which the wetting phase reconnects and displaces
    // the non-wetting phase back out of an already-invaded throat (imbibition).
    // [[codeblock]]
    auto getPcEntry = [&](const std::size_t eIdx)
    {
        const Scalar throatRadius = gridGeometry->throatInscribedRadius(eIdx);
        const auto shapeFactor = gridGeometry->throatShapeFactor(eIdx);
        return PoreNetwork::ThresholdCapillaryPressures::pcEntry(surfaceTension,
                                                                 contactAngle,
                                                                 throatRadius,
                                                                 shapeFactor);
    };

    auto getPcSnapOff = [&](const std::size_t eIdx)
    {
        const Scalar throatRadius = gridGeometry->throatInscribedRadius(eIdx);
        const auto shape = gridGeometry->throatCrossSectionShape(eIdx);
        return PoreNetwork::ThresholdCapillaryPressures::pcSnapoff(surfaceTension,
                                                                    contactAngle,
                                                                    throatRadius,
                                                                    shape);
    };
    // [[/codeblock]]

    // ### Allocate simulation data and set up VTK output
    // `elementIsInvaded`/`elementIsTrapped` are the per-throat invasion/trapping state
    // tracked by `TwoPStatic`; `pc`/`sw` are the per-pore local capillary pressure and
    // wetting-phase saturation, which are what is averaged into the pc-S curve.
    // [[codeblock]]
    std::vector<bool> elementIsInvaded(leafGridView.size(0), false);
    std::vector<bool> elementIsTrapped(leafGridView.size(0), false);
    std::vector<bool> elementWasTrapped(leafGridView.size(0), false);
    std::vector<Scalar> pcEntry(leafGridView.size(0));
    std::vector<Scalar> pcSnapOff(leafGridView.size(0));
    std::vector<int> throatLabel(leafGridView.size(0));
    std::vector<Scalar> poreVolume(leafGridView.size(1));
    std::vector<Scalar> pc(leafGridView.size(1), 0.0);
    std::vector<Scalar> sw(leafGridView.size(1), 0.0);
    std::vector<int> poreLabel(leafGridView.size(1));

    static const auto name = getParam<std::string>("Problem.Name");
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    auto writer = std::make_shared<Dune::VTKWriter<GridView>>(leafGridView);
    Dune::VTKSequenceWriter<GridView> sequenceWriter(writer, name);
    sequenceWriter.addCellData(pcEntry, "pcEntry");
    sequenceWriter.addCellData(pcSnapOff, "pcSnapOff");
    sequenceWriter.addCellData(elementIsInvaded, "invaded");
    sequenceWriter.addCellData(elementIsTrapped, "trapped");
    sequenceWriter.addCellData(throatLabel, "throatLabel");
    sequenceWriter.addVertexData(poreLabel, "poreLabel");
    sequenceWriter.addVertexData(pc, "pc");
    sequenceWriter.addVertexData(sw, "sw");
    sequenceWriter.addVertexData(poreVolume, "poreVolume");

    for (const auto& element : elements(leafGridView))
    {
        const auto eIdx = leafGridView.indexSet().index(element);
        pcEntry[eIdx] = getPcEntry(eIdx);
        pcSnapOff[eIdx] = getPcSnapOff(eIdx);
        throatLabel[eIdx] = gridGeometry->throatLabel(eIdx);

        for (int i = 0; i < 2; ++i)
        {
            const auto dofIdx = leafGridView.indexSet().subIndex(element, i, 1 /*dofCodim*/);
            poreLabel[dofIdx] = gridGeometry->poreLabel(dofIdx);
            poreVolume[dofIdx] = gridGeometry->poreVolume(dofIdx);
        }
    }
    // [[/codeblock]]

    // ### Set up the quasi-static two-phase pore-network model
    // `TwoPStatic` performs the invasion-percolation sweep (with snap-off) for a given
    // global capillary pressure; see [`staticpnm.hh`](../../dumux/porenetwork/2p/static/staticpnm.hh).
    // [[codeblock]]
    const Scalar pcRange = finalPc - initialPc;
    const Scalar deltaPc = pcRange / numSteps;
    Scalar pcGlobal = initialPc;

    PoreNetwork::TwoPStatic<GridGeometry, Scalar> staticModel(*gridGeometry,
                                                               pcEntry,
                                                               pcSnapOff,
                                                               throatLabel,
                                                               inletThroatLabel,
                                                               outletThroatLabel,
                                                               allowDraingeOfOutlet);

    std::ofstream logfile;
    const auto logfileName = name + "_pc-s-curve.txt";
    logfile.open(logfileName);
    Scalar averageSaturation = 0;

    // Sw/pc history of every step, kept separately from the logfile so the drainage leg
    // (the first numSteps+1 entries) and the imbibition leg (the remaining entries,
    // sharing the apex point) can be plotted with distinct styles further down.
    std::vector<Scalar> swHistory;
    std::vector<Scalar> pcHistory;
    swHistory.reserve(2*numSteps + 1);
    pcHistory.reserve(2*numSteps + 1);
    const Scalar totalPoreVolume = [&]()
    {
        Scalar result = 0.0;
        for (std::size_t i = 0; i < poreVolume.size(); ++i)
            if (poreLabel[i] != inletPoreLabel && (poreLabel[i] != outletPoreLabel || allowDraingeOfOutlet))
                result += poreVolume[i];
        return result;
    }();

    std::cout << "total pore volume is " << totalPoreVolume << std::endl;
    // [[/codeblock]]

    // ### The drainage-imbibition loop
    // At each step we (a) update the throat invasion state for the current global
    // capillary pressure, (b) update each pore's local capillary pressure/saturation --
    // gated on the *previous* step's trapped state, so that a pore which becomes newly
    // trapped this very step still receives its last legitimate update -- and (c) only
    // then recompute the trapped state for the next step. Swapping (b) and (c) would
    // silently discard that last saturation increment; see the class documentation of
    // `TwoPStatic::updateTrappedState` for the full reasoning.
    // [[codeblock]]
    for (int step = 0; step < 2*numSteps + 1; ++step)
    {
        std::cout << "Step " << step << ": Applying global pc of " << pcGlobal << " --> ";
        staticModel.updateInvasionState(elementIsInvaded, elementIsTrapped, pcGlobal);

        averageSaturation = 0;
        auto fvGeometry = localView(*gridGeometry);
        std::fill(sw.begin(), sw.end(), 1.0);

        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);
            elementWasTrapped[eIdx] = elementIsTrapped[eIdx];

            // pores connected to the invading phase receive the global capillary
            // pressure; trapped throats are excluded since their local capillary
            // pressure is decoupled from pcGlobal once trapped
            if (elementIsInvaded[eIdx] && !elementWasTrapped[eIdx])
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

            // get the saturation of each pore by inverting its local pc-s curve
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdx = scv.dofIndex();
                if (poreLabel[dofIdx] == inletPoreLabel || (poreLabel[dofIdx] == outletPoreLabel && !allowDraingeOfOutlet))
                    continue;

                if (pc[dofIdx] > 0.0)
                {
                    const Scalar poreRadius = gridGeometry->poreInscribedRadius(dofIdx);
                    const auto poreShape = gridGeometry->poreGeometry(dofIdx);
                    sw[dofIdx] = poreSaturation(poreShape, poreRadius, surfaceTension, pc[dofIdx]);
                }
                else
                    sw[dofIdx] = 1.0;

                const Scalar partialPoreVolume = scv.volume();
                averageSaturation += partialPoreVolume*sw[dofIdx];
            }
        }
        averageSaturation /= totalPoreVolume;

        // recompute which throats are trapped (disconnected from the inlet) for the
        // next step, now that this step's saturation update has already consumed the
        // previous trapped state
        staticModel.updateTrappedState(elementIsTrapped, elementIsInvaded);

        logfile << averageSaturation << " " << pcGlobal << std::endl;
        swHistory.push_back(averageSaturation);
        pcHistory.push_back(pcGlobal);

        if (step < numSteps)
            pcGlobal += deltaPc;
        else
            pcGlobal -= deltaPc;

        std::cout << staticModel.numThroatsInvaded() << " of " << elementIsInvaded.size()
                   << " throats invaded; S_avg " << averageSaturation << std::endl;
        sequenceWriter.write(step);
    }
    // [[/codeblock]]

    // ### Plot the resulting pc-S curve
    // The drainage leg (the first `numSteps+1` entries of the history, increasing `pcGlobal`)
    // and the imbibition leg (the remaining entries, decreasing `pcGlobal`, sharing the apex
    // point with the drainage leg) are plotted as two separate, distinctly styled and labeled
    // curves, so the hysteresis loop is visually unambiguous rather than looking like a single
    // folded-over line.
    // [[codeblock]]
#ifdef DUMUX_HAVE_GNUPLOT
    if (getParam<bool>("Problem.PlotPcS"))
    {
        const auto apex = static_cast<std::size_t>(numSteps);
        const std::vector<Scalar> swDrainage(swHistory.begin(), swHistory.begin() + apex + 1);
        const std::vector<Scalar> pcDrainage(pcHistory.begin(), pcHistory.begin() + apex + 1);
        const std::vector<Scalar> swImbibition(swHistory.begin() + apex, swHistory.end());
        const std::vector<Scalar> pcImbibition(pcHistory.begin() + apex, pcHistory.end());

        Dumux::GnuplotInterface<Scalar> gnuplot(true);
        gnuplot.setOpenPlotWindow(true);
        gnuplot.addDataSetToPlot(swDrainage, pcDrainage, name + "_drainage.dat",
                                 "with linespoints lt 1 lc rgb '#1f77b4' pt 7 ps 0.6 title 'Drainage'");
        gnuplot.addDataSetToPlot(swImbibition, pcImbibition, name + "_imbibition.dat",
                                 "with linespoints lt 2 lc rgb '#d62728' pt 5 ps 0.6 title 'Imbibition'");
        gnuplot.setXlabel("S_w [-]");
        gnuplot.setYlabel("p_c [Pa]");
        gnuplot.setXRange(0, 1.01);
        gnuplot.plot("plot");
    }
#endif
    // [[/codeblock]]

    // print dumux end message
    // [[codeblock]]
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
// [[/codeblock]]
// [[/details]]
// [[/content]]
