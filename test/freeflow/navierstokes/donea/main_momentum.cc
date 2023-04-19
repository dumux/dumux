// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model (Donea 2003, \cite Donea2003).
 */

#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/evalsolution.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/vtk/intersectionwriter.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>

#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "properties_momentum.hh"

namespace Dumux {

template<class Error>
void writeError_(std::ofstream& logFile, const Error& error)
{
    for (const auto& e : error)
        logFile << Fmt::format(", {:.5e}", e);
}

template<class Problem, class GridVariables, class SolutionVector>
void printErrors(std::shared_ptr<Problem> problem,
                 const GridVariables& gridVariables,
                 const SolutionVector& x)
{
    using GridGeometry = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    static constexpr int dim = GridGeometry::GridView::dimension;
    const bool printErrors = getParam<bool>("Problem.PrintErrors", false);

    if (printErrors)
    {
        NavierStokesTest::ErrorsSubProblem errors(problem, x);

        std::ofstream logFile(problem->name() + ".csv", std::ios::app);
        auto totalVolume = errors.totalVolume();
        auto numDofs = errors.numDofs();
        // For the staggered scheme, the control volumes are overlapping
        if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcstaggered)
            totalVolume /= dim;

        logFile << Fmt::format("{:.5e}", errors.time()) << ", ";
        logFile << numDofs << ", ";
        logFile << std::pow(totalVolume / numDofs, 1.0/dim);
        const auto& componentErrors = errors.l2Absolute();
        // Calculate L2-error for velocity field
        Dune::FieldVector<double, 1> velError(0.0);
        velError[0] = std::sqrt(componentErrors * componentErrors);

        writeError_(logFile, velError);
        logFile << "\n";
    }
}

template<class GridGeometry, class GridVariables, class SolutionVector>
void updateVelocities(
    std::vector<Dune::FieldVector<double, 2>>& velocity,
    std::vector<Dune::FieldVector<double, 2>>& faceVelocity,
    const GridGeometry& gridGeometry,
    const GridVariables& gridVariables,
    const SolutionVector& x
){
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVariables.curGridVolVars());
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);
        const auto eIdx = gridGeometry.elementMapper().index(element);

        if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcstaggered)
        {
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& vars = elemVolVars[scv];
                velocity[eIdx][scv.dofAxis()] += 0.5*vars.velocity();
                faceVelocity[scv.dofIndex()][scv.dofAxis()] = vars.velocity();
            }
        }
        else if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcdiamond)
        {
            const auto elemGeo = element.geometry();
            const auto elemSol = elementSolution(element, x, gridGeometry);
            velocity[eIdx] = evalSolution(element, elemGeo, gridGeometry, elemSol, elemGeo.center());
            for (const auto& scv : scvs(fvGeometry))
                faceVelocity[scv.dofIndex()] = elemVolVars[scv].velocity();
        }
        else if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::pq1bubble)
        {
            const auto elemGeo = element.geometry();
            const auto elemSol = elementSolution(element, x, gridGeometry);
            velocity[eIdx] = evalSolution(element, elemGeo, gridGeometry, elemSol, elemGeo.center());
        }
        else
            DUNE_THROW(Dune::Exception, "Unknown discretization type: " << GridGeometry::discMethod);
    }
}

template<class GridGeometry>
void updateRank(
    std::vector<int>& rank,
    const GridGeometry& gridGeometry
){
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        const auto eIdxGlobal = gridGeometry.elementMapper().index(element);
        rank[eIdxGlobal] = gridGeometry.gridView().comm().rank();
    }
}

} // end namespace Dumux


int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::DoneaTestMomentum;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    using GridManager = Dumux::GridManager<GetPropType<TypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // solve Stokes problem on this grid
    ////////////////////////////////////////////////////////////

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    std::vector<Dune::FieldVector<double, 2>> velocity(gridGeometry->gridView().size(0));
    std::vector<Dune::FieldVector<double, 2>> faceVelocity(x.size());
    updateVelocities(velocity, faceVelocity, *gridGeometry, *gridVariables, x);

    std::vector<int> rank(gridGeometry->gridView().size(0));
    updateRank(rank, *gridGeometry);

    std::vector<std::size_t> dofIdx(x.size());
    std::iota(dofIdx.begin(), dofIdx.end(), 0);

    std::string baseName = problem->name();
    std::string discSuffix = std::string("_") + GridGeometry::DiscretizationMethod::name();
    std::string rankSuffix = std::string("_") + std::to_string(gridGeometry->gridView().comm().rank());
    Dune::VTKWriter<typename GridGeometry::GridView> writer(gridGeometry->gridView());
    using Field = Vtk::template Field<typename GridGeometry::GridView>;
    writer.addCellData(Field(
        gridGeometry->gridView(), gridGeometry->elementMapper(), velocity,
        "velocity", /*numComp*/2, /*codim*/0
    ).get());
    writer.addCellData(rank, "rank");
    writer.write(baseName + discSuffix + "_0");

    ConformingIntersectionWriter faceVtk(gridGeometry->gridView());
    // face quantities have no special significance for the PQ1Bubble scheme
    if constexpr (GridGeometry::discMethod != DiscretizationMethods::pq1bubble)
    {
        faceVtk.addField(dofIdx, "dofIdx");
        faceVtk.addField(faceVelocity, "velocityVector");
        faceVtk.addField([&](const auto& is, const auto idx) {
            const auto& facet = is.inside().template subEntity <1> (is.indexInInside());
            return facet.partitionType();
        }, "partitionType");
        faceVtk.write(baseName + "_face" + discSuffix + rankSuffix + "_0", Dune::VTK::ascii);
    }

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<GridGeometry>,
                                                  LinearAlgebraTraitsFromAssembler<Assembler>>;
    const auto& dofMapper = LinearSolverTraits<GridGeometry>::dofMapper(*gridGeometry);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), dofMapper);

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    ////////////////////////////////////////////////////////////
    // write VTK output
    ////////////////////////////////////////////////////////////
    Dumux::updateVelocities(velocity, faceVelocity, *gridGeometry, *gridVariables, x);
    writer.write(baseName + discSuffix + "_1");

    if constexpr (GridGeometry::discMethod != DiscretizationMethods::pq1bubble)
        faceVtk.write(baseName + "_face" + discSuffix + rankSuffix + "_1", Dune::VTK::ascii);

    Dumux::printErrors(problem, *gridVariables, x);

    ////////////////////////////////////////////////////////////
    // finalize, print parameters and Dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
