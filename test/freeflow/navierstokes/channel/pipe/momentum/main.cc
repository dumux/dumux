// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Pipe flow Stokes test with known pressure
 */

#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <test/freeflow/navierstokes/errors.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    // get an instance of the MPI helper
    const auto& mpiHelper = Dune::MPIHelper::instance();

    Parameters::init(argc, argv);

    using GridManager = Dumux::GridManager<GetPropType<TypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

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

    const auto vtkOutputType = [&]{
        if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcdiamond)
            return Dune::VTK::nonconforming;
        else
            return Dune::VTK::conforming;
    }();

    using VTKOut = VtkOutputModule<GridVariables, SolutionVector>;
    VTKOut vtkWriter(*gridVariables, x, problem->name(), "", vtkOutputType);

    vtkWriter.addVolumeVariable([](const auto& v){ return v.velocity(); }, "velocity");
    std::vector<double> pressure(gridGeometry->gridView().size(0));
    for (const auto& e : elements(gridGeometry->gridView()))
    {
        const auto center = e.geometry().center();
        pressure[gridGeometry->elementMapper().index(e)] = problem->analyticalPressure(center);
    }
    vtkWriter.addField(pressure, "p");
    vtkWriter.write(0.0);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    vtkWriter.write(1.0);

    NavierStokesTest::ErrorsSubProblem<Problem> errors(problem, x);
    convergenceTestAppendErrorsMomentum(problem, errors);

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
}
