// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::CooksMembrane;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
    std::cout << "Number of Dirichlet constraints: " << problem->constraints().size() << "\n";

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addVolumeVariable([](const auto& v){
        return Dune::FieldVector<double, 2>{ v.displacement(0), v.displacement(1) };
    }, "u");
    vtkWriter.write(0.0);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // PQ2 with nearly-incompressible materials produces an ill-conditioned stiffness
    // matrix (condition number ∝ λ/µ). AMG struggles with this; use UMFPack instead.
    using DiscMethod = typename GridGeometry::DiscretizationMethod;
    using LinearSolver = std::conditional_t<
        std::is_same_v<DiscMethod, DiscretizationMethods::PQ2>,
        UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>,
        AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>
    >;
    auto linearSolver = [&]() -> std::shared_ptr<LinearSolver> {
        if constexpr (std::is_same_v<DiscMethod, DiscretizationMethods::PQ2>)
            return std::make_shared<LinearSolver>();
        else
            return std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
    }();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // Debug: check residual and Jacobian before Newton
    {
        assembler->assembleJacobianAndResidual(x);
        const auto& res = assembler->residual();
        const auto& jac = assembler->jacobian();
        const int numEq = 2;
        std::cout << "\n[DEBUG] numDofs=" << x.size() << "\n";

        // Count DOF types (vertices vs edge midpoints) using the DOF mapper
        const auto& gg = *gridGeometry;
        const auto& dofMapper = gg.dofMapper();
        int nVertexDofs = 0, nEdgeDofs = 0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            // Check if dof i is a vertex or edge midpoint dof
            // by checking if it appears in the vertex mapper
            // Use the grid's vertex and edge index sets
        }

        // Check non-Dirichlet (free) rows: if row norm is ~1 it's Dirichlet
        // otherwise it's a free DOF
        int freeDofCount = 0;
        int dirichletDofCount = 0;
        double minFreeRowNorm = 1e300;
        std::size_t minFreeRow = 0; int minFreeEq = 0;

        for (std::size_t i = 0; i < x.size(); ++i) {
            for (int eq = 0; eq < numEq; ++eq) {
                // check if Dirichlet: only diagonal [eq][eq] = 1, rest = 0
                bool isDirichlet = true;
                const auto& row = jac[i];
                for (auto col = row.begin(); col != row.end(); ++col) {
                    for (int c = 0; c < numEq; ++c) {
                        double val = (*col)[eq][c];
                        bool onDiag = (col.index() == i && c == eq);
                        double expected = onDiag ? 1.0 : 0.0;
                        if (std::abs(val - expected) > 1e-10) { isDirichlet = false; break; }
                    }
                    if (!isDirichlet) break;
                }
                if (isDirichlet) { dirichletDofCount++; continue; }
                freeDofCount++;

                double rowNorm = 0.0;
                for (auto col = row.begin(); col != row.end(); ++col)
                    for (int c = 0; c < numEq; ++c) rowNorm += (*col)[eq][c] * (*col)[eq][c];
                rowNorm = std::sqrt(rowNorm);
                if (rowNorm < minFreeRowNorm) { minFreeRowNorm = rowNorm; minFreeRow = i; minFreeEq = eq; }
            }
        }

        std::cout << "[DEBUG] Dirichlet rows: " << dirichletDofCount
                  << ", free rows: " << freeDofCount << "\n";
        std::cout << "[DEBUG] min free row norm: " << minFreeRowNorm
                  << " at dof=" << minFreeRow << " eq=" << minFreeEq << "\n";
        std::cout << "[DEBUG] row at min free dof:\n";
        for (auto col = jac[minFreeRow].begin(); col != jac[minFreeRow].end(); ++col)
            std::cout << "  col=" << col.index() << " block=["
                      << (*col)[minFreeEq][0] << "," << (*col)[minFreeEq][1] << "]\n";

        double residNorm = 0.0;
        for (std::size_t i = 0; i < res.size(); ++i) residNorm += res[i].two_norm2();
        std::cout << "[DEBUG] residual norm: " << std::sqrt(residNorm) << "\n\n";
    }

    nonLinearSolver.solve(x);
    gridVariables->update(x);
    vtkWriter.write(1.0);

    // report tip displacement at the top-right corner (48, 60)
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    const GlobalPosition tipPoint{ 48.0, 60.0 };
    const auto tipElements = intersectingEntities(tipPoint, gridGeometry->boundingBoxTree());
    if (!tipElements.empty())
    {
        const auto element = gridGeometry->element(tipElements[0]);
        const auto elemSol = elementSolution(element, x, *gridGeometry);
        const auto u = evalSolution(element, element.geometry(), *gridGeometry, elemSol, tipPoint);
        std::cout << "Tip displacement at (48, 60):"
                  << "  u_x = " << u[0]
                  << "  u_y = " << u[1] << "\n";
    }

    return 0;
}
