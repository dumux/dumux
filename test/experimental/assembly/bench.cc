// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief What analytic differentiation is worth, by scheme
 *
 * Numeric differentiation re-evaluates the element residual once per local dof and equation, so
 * its cost grows with the size of the local basis. An analytic Jacobian is one extra pass over the
 * element. The interesting question is therefore not the ratio on any single scheme but how it
 * *grows*, which is what the hybrid schemes — 6, 9 and 16 local dofs in 2D against box's 4 — are
 * here to show.
 *
 * The residual column is the same code in both configurations; it is the floor neither method can
 * go below, and it is what makes the ratio less than the deflection count.
 */
#include <config.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/localview.hh>

#include "properties.hh"

namespace {

template<class Function>
double best(Function&& run, int repeats = 8)
{
    run();
    double fastest = std::numeric_limits<double>::infinity();
    for (int i = 0; i < repeats; ++i)
    {
        const auto start = std::chrono::steady_clock::now();
        run();
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        fastest = std::min(fastest, elapsed.count());
    }
    return fastest;
}

template<class TypeTag>
void bench(const std::string& name, int cells)
{
    using Grid = Dumux::GetPropType<TypeTag, Dumux::Properties::Grid>;
    using GridGeometry = Dumux::GetPropType<TypeTag, Dumux::Properties::GridGeometry>;
    using Problem = Dumux::GetPropType<TypeTag, Dumux::Properties::Problem>;
    using GridVariables = Dumux::GetPropType<TypeTag, Dumux::Properties::GridVariables>;
    using SolutionVector = Dumux::GetPropType<TypeTag, Dumux::Properties::SolutionVector>;

    Dune::FieldVector<double, 2> upperRight(1.0);
    std::array<int, 2> numCells{cells, cells};
    Grid grid(upperRight, numCells);

    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    SolutionVector x(gridGeometry->numDofs());
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i] = 0.3 + 0.4*std::sin(1.7*double(i)) + 0.2*std::cos(0.9*double(i));

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // the number of local dofs is what sets how often numeric differentiation repeats itself
    const auto fvGeometry = localView(*gridGeometry).bindElement(*elements(grid.leafGridView()).begin());
    const auto numLocalDofs = Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry);

    Dumux::Experimental::Assembler<TypeTag, Dumux::DiffMethod::analytic> analytic(
        problem, gridGeometry, gridVariables);
    Dumux::Experimental::Assembler<TypeTag, Dumux::DiffMethod::numeric> numeric(
        problem, gridGeometry, gridVariables);
    analytic.setLinearSystem();
    numeric.setLinearSystem();

    const auto tAnalytic = best([&]{ analytic.assembleJacobianAndResidual(x); });
    const auto tNumeric = best([&]{ numeric.assembleJacobianAndResidual(x); });
    const auto tResidual = best([&]{ analytic.assembleResidual(x); });

    std::cout << std::left << std::setw(19) << name << std::right
              << std::setw(8) << gridGeometry->numDofs()
              << std::setw(7) << numLocalDofs
              << std::setw(11) << std::fixed << std::setprecision(1) << tAnalytic*1e3
              << std::setw(10) << tNumeric*1e3
              << std::setw(10) << tResidual*1e3
              << std::setw(9) << std::setprecision(2) << tNumeric/tAnalytic << "x"
              << std::endl;

    // machine-readable too: the run-to-run spread on a loaded machine is larger than any effect
    // worth reporting, so quote the minimum over repeated *processes*, not over one process
    std::cout << "DATA;" << name << ";" << gridGeometry->numDofs() << ";" << numLocalDofs
              << ";" << tAnalytic*1e3 << ";" << tNumeric*1e3 << ";" << tResidual*1e3
              << std::endl;
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init([](Dune::ParameterTree& params){
        params["Assembly.NumericDifferenceMethod"] = "0";
        params["Assembly.NumericDifference.BaseEpsilon"] = "1e-6";
        // one thread, so the comparison is of arithmetic rather than of scheduling
        params["Assembly.Multithreading"] = "false";
    });

    std::cout << "Analytic vs numeric differentiation, experimental CVFE assembler\n"
              << std::left << std::setw(19) << "scheme" << std::right
              << std::setw(8) << "dofs" << std::setw(7) << "local"
              << std::setw(11) << "analytic" << std::setw(10) << "numeric"
              << std::setw(10) << "residual" << std::setw(10) << "speedup" << std::endl;

    bench<Dumux::Properties::TTag::NonlinearDiffusionPQ1Bubble>("pq1bubble", 64);
    bench<Dumux::Properties::TTag::NonlinearDiffusionPQ1BubbleHybrid>("pq1bubble hybrid", 64);
    bench<Dumux::Properties::TTag::NonlinearDiffusionPQ2Hybrid>("pq2 hybrid", 48);
    bench<Dumux::Properties::TTag::NonlinearDiffusionPQ3Hybrid>("pq3 hybrid", 32);

    return 0;
}
