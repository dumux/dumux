// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Drive the parallel FoamGrid through the standard Dumux::GridManager: GridManager reads the
 *        network and (opt-in via Grid.ParallelPartitioning) distributes it spatially; we then solve
 *        a 1p cc-tpfa problem on the distributed grid and check the partition and the solve.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

#include "properties.hh"

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

//! minimal data handle to communicate one double per element (owner -> overlap copies)
template<class GridView>
class ElementScalarHandle
: public Dune::CommDataHandleIF<ElementScalarHandle<GridView>, double>
{
public:
    ElementScalarHandle(const GridView& gv, std::vector<double>& data) : gv_(gv), data_(data) {}
    bool contains(int, int codim) const { return codim == 0; }
    bool fixedSize(int, int) const { return true; }
    template<class E> std::size_t size(const E&) const { return 1; }
    template<class Buf, class E> void gather(Buf& b, const E& e) const { b.write(data_[gv_.indexSet().index(e)]); }
    template<class Buf, class E> void scatter(Buf& b, const E& e, std::size_t) { double v; b.read(v); data_[gv_.indexSet().index(e)] = v; }
private:
    const GridView& gv_;
    std::vector<double>& data_;
};

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::OnePNetworkCCTpfa;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;

    // the GridManager reads the network and (opt-in) distributes it over the ranks
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    // the grid view exposes the distributed communicator after the GridManager has load-balanced
    const auto& comm = leafGridView.comm();
    int failures = 0;

    // the network.dgf has 1735 elements; each must be interior on exactly one rank
    int nInterior = 0;
    for (const auto& e : elements(leafGridView))
        if (e.partitionType() == Dune::InteriorEntity) ++nInterior;
    const int total = comm.sum(nInterior);
    if (total != 1735) { if (comm.rank() == 0) std::cerr << "FAIL: interior sum " << total << " != 1735\n"; ++failures; }
    if (comm.size() > 1 && leafGridView.overlapSize(0) == 0) { if (comm.rank() == 0) std::cerr << "FAIL: no overlap\n"; ++failures; }

    // assemble and solve a 1p problem on the distributed grid
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;

    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    auto problem = std::make_shared<Problem>(gridGeometry);
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
    NewtonSolver<Assembler, LinearSolver> newton(assembler, linearSolver);
    newton.solve(x);

    // VTK output: in parallel this writes one .vtu per rank plus a .pvtu, exercising the grid's
    // partition types (only interior/border entities are written) and global id set
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, "test_1p_network1d3d_gridmanager");
    IOFields::initOutputModule(vtkWriter);
    std::vector<int> rank(leafGridView.size(0), comm.rank()); // make the decomposition visible
    vtkWriter.addField(rank, "rank");
    vtkWriter.write(0.0);

    // the owned solution must be finite
    using std::isfinite;
    int nOwnedFinite = 0;
    for (const auto& e : elements(leafGridView, Dune::Partitions::interior))
        if (isfinite(x[gridGeometry->elementMapper().index(e)][0])) ++nOwnedFinite;
    if (comm.sum(nOwnedFinite) != 1735) { if (comm.rank() == 0) std::cerr << "FAIL: non-finite solution\n"; ++failures; }

    // the DGF per-element parameter (radius) must be migrated onto the distributed local elements
    // (incl. the overlap). Correctness: the radius an overlap element reports via the standard DuMux
    // GridData parameter migration (keyed by the now-persistent localIdSet id) must equal the
    // owner's radius, which we independently transport into the overlap copies with communicate().
    {
        auto gridData = gridManager.getGridData();
        std::vector<double> ownerRadius(leafGridView.size(0), -1.0);
        for (const auto& e : elements(leafGridView, Dune::Partitions::interior))
            ownerRadius[gridGeometry->elementMapper().index(e)] = gridData->parameters(e)[0];
        if (comm.size() > 1)
        {
            ElementScalarHandle<std::decay_t<decltype(leafGridView)>> handle(leafGridView, ownerRadius);
            leafGridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
        }

        int nBadRadius = 0;
        for (const auto& e : elements(leafGridView))
        {
            const double migrated = gridData->parameters(e)[0];
            const double viaComm = ownerRadius[gridGeometry->elementMapper().index(e)];
            if (!(migrated > 0.0) || !isfinite(migrated) || std::abs(migrated - viaComm) > 1e-12*std::abs(viaComm))
                ++nBadRadius;
        }
        if (comm.sum(nBadRadius) != 0) { if (comm.rank() == 0) std::cerr << "FAIL: DGF radius not migrated correctly\n"; ++failures; }
    }

    const int totalFailures = comm.sum(failures);
    if (comm.rank() == 0)
        std::cout << (totalFailures == 0
            ? "GridManager distributed + solved the 1p network: PASS\n"
            : "GridManager parallel 1p: FAIL\n");
    return totalFailures > 0 ? 1 : 0;
}
