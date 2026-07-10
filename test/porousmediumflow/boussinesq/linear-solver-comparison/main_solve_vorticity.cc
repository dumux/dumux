// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Matrix-solving stage for the vector-potential (vorticity) Boussinesq
 *        linear-solver benchmark.
 *
 * Reads the Jacobian/residual pairs written by main_generate_vorticity.cc
 * (one per time in MatrixExport.Times) from MatrixMarket files and
 * benchmarks each candidate ISTL solver against them (see solvercompare.hh).
 * No grid, problem or transient simulation is set up here -- solving is
 * fully decoupled from how the matrices were produced.
 */
#include <config.h>

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/matrixmarket.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/assembly/fvassembler.hh>

#include "properties_vorticity.hh"
#include "matrixio.hh"
#include "solvercompare.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::BoussinesqOneSidedRB;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv, [](auto&){}, "params.input");

    using Scalar    = GetPropType<TypeTag, Properties::Scalar>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using Matrix    = typename Assembler::JacobianMatrix;
    using Vector    = typename Assembler::ResidualType;

    const auto problemName = getParam<std::string>("Problem.Name");
    const auto exportTimes = getParam<std::vector<Scalar>>("MatrixExport.Times");

    std::ofstream csv(problemName + "_linsolver_comparison.csv");
    writeSolverComparisonCsvHeader(csv);

    for (const auto t : exportTimes)
    {
        Matrix A;
        Vector b;
        Dune::loadMatrixMarket(A, matrixFileName(problemName, t));
        Dune::loadMatrixMarket(b, rhsFileName(problemName, t));

        const std::string state = "t=" + std::to_string(static_cast<long long>(std::llround(t)));
        const auto results = benchmarkAllSolvers<Assembler>(state, A, b);
        printSolverComparison("Vorticity formulation -- " + state, results);
        writeSolverComparisonCsv(csv, "vorticity", results);
    }

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
