// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Convergence test with analytic solution
 */
#include <config.h>
#include <iostream>
#include <iomanip>
#include <type_traits>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/gridwriter.hh>
#include <dumux/io/cvfegridfunction.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/discretization/concepts.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/assembler.hh>

#include <dumux/common/metadata.hh>

#include "properties.hh"

template<class Problem, class SolutionVector>
void printL2Error(const Problem& problem, const SolutionVector& x)
{
    using namespace Dumux;
    using Scalar = double;

    Scalar l2error = 0.0;
    auto fvGeometry = localView(problem.gridGeometry());
    for (const auto& element : elements(problem.gridGeometry().gridView()))
    {
        fvGeometry.bindElement(element);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto dofIdx = scv.dofIndex();
            const Scalar delta = x[dofIdx] - problem.analyticalSolution(scv.dofPosition())[2/*pressureIdx*/];
            l2error += scv.volume()*(delta*delta);
        }
    }
    using std::sqrt;
    l2error = sqrt(l2error);

    const auto numDofs = problem.gridGeometry().numDofs();
    std::ostream tmp(std::cout.rdbuf());
    tmp << std::setprecision(8) << "** L2 error (abs) for "
            << std::setw(6) << numDofs << " dofs "
            << std::scientific
            << "L2 error = " << l2error
            << std::endl;

    Dumux::MetaData::Collector collector;
    if (Dumux::MetaData::jsonFileExists(problem.name()))
        Dumux::MetaData::readJsonFile(collector, problem.name());
    collector["numDofs"].push_back(numDofs);
    collector["L2errors"].push_back(l2error);
    Dumux::MetaData::writeJsonFile(collector, problem.name());
}

//! whether the variables are defined per local dof rather than per sub-control volume
template<class GridVariables>
constexpr bool usesGeneralGridVariables()
{ return Dumux::Concept::GridVariables<GridVariables> && !Dumux::Concept::FVGridVariables<GridVariables>; }

/*!
 * \brief Write the solution to vtk
 * \note The variables defined per local dof are not known to the output module, so the
 *       solution is written with the grid writer at the interpolation order of the scheme.
 */
template<class TypeTag, class Problem, class GridGeometry, class GridVariables, class SolutionVector>
void writeVtkOutput(const Problem& problem,
                    const GridGeometry& gridGeometry,
                    GridVariables& gridVariables,
                    const SolutionVector& x)
{
    using namespace Dumux;

    if constexpr (usesGeneralGridVariables<GridVariables>())
    {
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::pq3)
        {
            IO::GridWriter writer{IO::Format::vtu, gridGeometry.gridView(), IO::order<3>};
            writer.setPointField("p", x);
            writer.write(problem.name());
        }
        else
        {
            IO::GridWriter writer{IO::Format::vtu, gridGeometry.gridView(), IO::order<2>};
            writer.setPointField("p", IO::cvfeGridFunction(gridGeometry, x));
            writer.write(problem.name());
        }
    }
    else
    {
        VtkOutputModule<GridVariables, SolutionVector> vtkWriter(gridVariables, x, problem.name());
        using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
        vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(gridVariables));
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        IOFields::initOutputModule(vtkWriter); // Add model specific output fields
        vtkWriter.write(1.0);
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

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
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // create assembler & linear solver
    using Assembler = std::conditional_t<usesGeneralGridVariables<GridVariables>(),
                                         Dumux::Experimental::Assembler<TypeTag, DiffMethod::numeric>,
                                         FVAssembler<TypeTag, DiffMethod::numeric>>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // linearize & solve
    nonLinearSolver.solve(x);

    writeVtkOutput<TypeTag>(*problem, *gridGeometry, *gridVariables, x);

    printL2Error(*problem, x);

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

}
