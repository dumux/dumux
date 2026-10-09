// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Cook's membrane with linear elasticity: displacement field and vertical displacement
 *        of the upper right corner
 */
#include <config.h>

#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#if DUMUX_HAVE_GRIDFORMAT
#include <dumux/io/gridwriter.hh>
#include <dumux/io/cvfegridfunction.hh>
#endif

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::CooksMembrane;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonlinearSolver(assembler, linearSolver);
    nonlinearSolver.solve(x);

    // the upper right corner is a vertex, where the discrete displacement is the value of its dof
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    const GlobalPosition tip{48.0, 60.0};
    auto elemDisc = localView(*gridGeometry);
    std::optional<std::size_t> tipDof;
    for (const auto& element : elements(leafGridView))
    {
        elemDisc.bind(element);
        for (const auto& localDof : localDofs(elemDisc))
            if ((ipData(elemDisc, localDof).global() - tip).two_norm() < 1e-8)
                tipDof = localDof.dofIndex();
        if (tipDof)
            break;
    }
    if (!tipDof)
        DUNE_THROW(Dune::InvalidStateException, "No degree of freedom at the tip (48, 60)");

    std::cout << "Number of dofs: " << gridGeometry->numDofs() << "\n"
              << std::setprecision(10) << "Tip displacement at (48, 60): u_x = " << x[*tipDof][0]
              << ", u_y = " << x[*tipDof][1] << std::endl;

#if DUMUX_HAVE_GRIDFORMAT
    using DiscretizationMethod = typename GridGeometry::DiscretizationMethod;
    if constexpr (std::is_same_v<DiscretizationMethod, DiscretizationMethods::Box>)
    {
        IO::GridWriter writer{IO::Format::vtu, leafGridView, IO::order<1>};
        writer.setPointField("u", IO::cvfeGridFunction(*gridGeometry, x));
        writer.write(problem->name());
    }
    else
    {
        IO::GridWriter writer{IO::Format::vtu, leafGridView, IO::order<2>};
        writer.setPointField("u", IO::cvfeGridFunction(*gridGeometry, x));
        writer.write(problem->name());
    }
#endif

    return 0;
}
