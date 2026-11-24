// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the analytic differentiation mode of the CVFE local assembler
 *
 * The same nonlinear diffusion problem is assembled twice, differing only in how the Jacobian is
 * obtained, and the two matrices are compared entry by entry. `DiffMethod::analytic` routes
 * through the local residual's `addElementJacobian`, `DiffMethod::numeric` through deflection of
 * every local dof, so agreement is a statement about the hand-written derivative and about the
 * assembler that calls it.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dune/common/parametertree.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/linear/dunevectors.hh>

#include "properties.hh"

namespace {

//! the largest relative difference between two assembled Jacobians
template<class Matrix>
double compare(const Matrix& analytic, const Matrix& numeric)
{
    double worst = 0.0;
    double scale = 0.0;
    for (auto row = numeric.begin(); row != numeric.end(); ++row)
        for (auto entry = row->begin(); entry != row->end(); ++entry)
        {
            const auto reference = (*entry)[0][0];
            const auto value = analytic[row.index()][entry.index()][0][0];
            scale = std::max(scale, std::abs(reference));
            worst = std::max(worst, std::abs(value - reference));
        }
    return worst/std::max(scale, 1e-30);
}

template<class TypeTag>
bool run(const std::string& name, int cells)
{
    using Grid = Dumux::GetPropType<TypeTag, Dumux::Properties::Grid>;
    using GridGeometry = Dumux::GetPropType<TypeTag, Dumux::Properties::GridGeometry>;
    using Problem = Dumux::GetPropType<TypeTag, Dumux::Properties::Problem>;
    using GridVariables = Dumux::GetPropType<TypeTag, Dumux::Properties::GridVariables>;
    using SolutionVector = Dumux::GetPropType<TypeTag, Dumux::Properties::SolutionVector>;

    Dune::FieldVector<double, 2> upperRight(1.0);
    std::array<int, 2> elements{cells, cells};
    Grid grid(upperRight, elements);

    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    // a linearisation point, not a solution: deliberately varying from dof to dof, so that
    // no two neighbouring shape-function contributions can cancel an error in the derivative
    SolutionVector x(gridGeometry->numDofs());
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] = 0.3 + 0.4*std::sin(1.7*double(i)) + 0.2*std::cos(0.9*double(i));

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Analytic = Dumux::Experimental::Assembler<TypeTag, Dumux::DiffMethod::analytic>;
    using Numeric = Dumux::Experimental::Assembler<TypeTag, Dumux::DiffMethod::numeric>;

    Analytic analytic(problem, gridGeometry, gridVariables);
    Numeric numeric(problem, gridGeometry, gridVariables);
    analytic.setLinearSystem();
    numeric.setLinearSystem();
    analytic.assembleJacobianAndResidual(x);
    numeric.assembleJacobianAndResidual(x);

    const auto jacobianError = compare(analytic.jacobian(), numeric.jacobian());

    auto residualDifference = analytic.residual();
    residualDifference -= numeric.residual();
    const auto residualError = residualDifference.infinity_norm()
                             / std::max(numeric.residual().infinity_norm(), 1e-30);

    // The bound on the Jacobian is a finite-difference bound, not a fudge. Sweeping the step
    // (pass it as argv[1]) traces the V-shape of central differences — pq3 hybrid, the sharpest
    // case, being
    //     eps  1e-4     1e-5     1e-6     1e-7     1e-8
    //     err  4.5e-9   4.2e-11  1.1e-10  1.2e-9   1.0e-8
    // truncation on the left, round-off on the right, minimum near 1e-11. A wrong analytic
    // derivative would floor at a constant instead: what converges here is the *numeric*
    // Jacobian, towards the analytic one. The residual is the same code in both and must agree
    // exactly.
    const bool passed = jacobianError < 1e-8 && residualError == 0.0;
    std::cout << name << ": " << gridGeometry->numDofs() << " dofs, "
              << "jacobian " << jacobianError << ", residual " << residualError
              << (passed ? "  PASS" : "  FAIL") << std::endl;
    return passed;
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    // runtime parameters, not defaults: `getParamFromGroup(group, key, fallback)` prefers its own
    // fallback over the default tree, so a default here would silently not apply
    Dumux::Parameters::init([argc, argv](Dune::ParameterTree& params){
        params["Assembly.NumericDifferenceMethod"] = "0";
        // overridable from the command line so the truncation claim below can be re-checked
        params["Assembly.NumericDifference.BaseEpsilon"] = argc > 1 ? argv[1] : "1e-6";
    });

    bool passed = true;
    passed &= run<Dumux::Properties::TTag::NonlinearDiffusionPQ1Bubble>("pq1bubble", 8);

    // the hybrid schemes are what an element-level Jacobian hook is really for: their
    // Galerkin-tested dofs have no sub-control volume, so no per-scv or per-scvf hook reaches them
    passed &= run<Dumux::Properties::TTag::NonlinearDiffusionPQ1BubbleHybrid>("pq1bubble hybrid", 8);
    passed &= run<Dumux::Properties::TTag::NonlinearDiffusionPQ2Hybrid>("pq2 hybrid", 6);
    passed &= run<Dumux::Properties::TTag::NonlinearDiffusionPQ3Hybrid>("pq3 hybrid", 4);

    std::cout << (passed ? "All tests passed" : "TEST FAILED") << std::endl;
    return passed ? 0 : 1;
}
