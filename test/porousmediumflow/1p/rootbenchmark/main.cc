// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Root benchmark case Schnepf et al 2020 M3.1 doi: 10.3389/fpls.2020.00316
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/assembly/fvassembler.hh>

#if DUMUX_HAVE_GNUPLOT
#include <dumux/io/gnuplotinterface.hh>
#endif
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addVolumeVariable([](const auto& v) { return v.pressure(); }, "p (Pa)");
    vtkWriter.addVolumeVariable([](const auto& v) { return 100*(v.pressure() - 1e5)/(1000*9.81); }, "psi (cm)");
    vtkWriter.addField(problem->analyticalPressure(), "p_exact");
    vtkWriter.addField(problem->analyticalHead(), "psi_exact");
    vtkWriter.write(0.0);

    // the system assembler
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the linear solver
    using LinearSolver = SSORCGIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the system solver
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;
    Solver solver(assembler, linearSolver);

    // solve & write output
    solver.solve(x);
    vtkWriter.write(1.0);

    // write more output for benchmark plot
    problem->outputL2Norm(x);

#if DUMUX_HAVE_GNUPLOT
    // plot with gnuplot
    bool enablePlot = getParam<bool>("Problem.EnablePlot", false);
    if (enablePlot)
    {
        GnuplotInterface<double> gnuplot;
        std::vector<double> psi(x.size()), z(x.size());
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            z[eIdx] = element.geometry().center()[2];
            psi[eIdx] = 100*(x[eIdx][0] - 1.0e5)/(1000*9.81);
        }
        gnuplot.addDataSetToPlot(psi, z, "psi.dat", "w l t 'DuMux'");
        gnuplot.addDataSetToPlot(problem->analyticalHead(), z, "psi_exact.dat", "w p t 'exact'");
        gnuplot.plot();
    }
#endif

    // print parameters (used and unused)
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

} // end main
