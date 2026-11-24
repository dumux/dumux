// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test that the experimental assembler also supports CCTpfa.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/fvassembler.hh>

#include "../../porousmediumflow/1p/incompressible/properties.hh"

namespace {

template<class MatrixA, class MatrixB>
double maxAbsDiff(const MatrixA& a, const MatrixB& b)
{
    double result = 0.0;
    for (auto row = a.begin(); row != a.end(); ++row)
        for (auto entry = row->begin(); entry != row->end(); ++entry)
        {
            const auto& otherBlock = b[row.index()][entry.index()];
            for (int i = 0; i < entry->N(); ++i)
                for (int j = 0; j < entry->M(); ++j)
                    result = std::max(result, std::abs((*entry)[i][j] - otherBlock[i][j]));
        }

    return result;
}

template<class VectorA, class VectorB>
double maxAbsDiffVector(const VectorA& a, const VectorB& b)
{
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a[i].size(); ++j)
            result = std::max(result, std::abs(a[i][j] - b[i][j]));

    return result;
}

template<Dumux::DiffMethod diffMethod>
void runComparison(const std::string& name)
{
    using TypeTag = Dumux::Properties::TTag::OnePIncompressibleTpfa;
    using Grid = Dumux::GetPropType<TypeTag, Dumux::Properties::Grid>;
    using GridGeometry = Dumux::GetPropType<TypeTag, Dumux::Properties::GridGeometry>;
    using Problem = Dumux::GetPropType<TypeTag, Dumux::Properties::Problem>;
    using GridVariables = Dumux::GetPropType<TypeTag, Dumux::Properties::GridVariables>;
    using SolutionVector = Dumux::GetPropType<TypeTag, Dumux::Properties::SolutionVector>;

    Dune::FieldVector<double, 2> upperRight(1.0);
    std::array<int, 2> cells{4, 4};
    Grid grid(upperRight, cells);

    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    SolutionVector x(gridGeometry->numDofs());
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        const auto dofIdx = gridGeometry->elementMapper().index(element);
        const auto center = element.geometry().center();
        x[dofIdx][0] = 1.0e5
                     + 2.0e4*std::sin(2.1*center[0])
                     - 1.0e4*std::cos(1.7*center[1]);
    }

    auto referenceGridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    auto experimentalGridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    referenceGridVariables->init(x);
    experimentalGridVariables->init(x);

    using ReferenceAssembler = Dumux::FVAssembler<TypeTag, diffMethod>;
    using ExperimentalAssembler = Dumux::Experimental::Assembler<TypeTag, diffMethod>;

    ReferenceAssembler reference(problem, gridGeometry, referenceGridVariables);
    ExperimentalAssembler experimental(problem, gridGeometry, experimentalGridVariables);
    reference.setLinearSystem();
    experimental.setLinearSystem();

    reference.assembleJacobianAndResidual(x);
    experimental.assembleJacobianAndResidual(x);

    const auto jacobianError = maxAbsDiff(reference.jacobian(), experimental.jacobian());
    const auto residualError = maxAbsDiffVector(reference.residual(), experimental.residual());

    std::cout << name << ": jacobian " << jacobianError
              << ", residual " << residualError << std::endl;

    if (jacobianError > 1e-20 || residualError > 1e-20)
        DUNE_THROW(Dune::InvalidStateException, "Experimental::Assembler differs from FVAssembler for CCTpfa");
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init([](Dune::ParameterTree& params)
    {
        params["Problem.Name"] = "experimental_assembler_tpfa";
        params["Problem.ExtrusionFactor"] = "1";
        params["Problem.CheckIsConstantVelocity"] = "false";
        params["SpatialParams.LensLowerLeft"] = "0.2 0.2";
        params["SpatialParams.LensUpperRight"] = "0.8 0.8";
        params["SpatialParams.Permeability"] = "1e-10";
        params["SpatialParams.PermeabilityLens"] = "1e-12";
        params["Assembly.NumericDifference.PriVarMagnitude"] = "1e5";
        params["Assembly.NumericDifference.BaseEpsilon"] = "0.1";
    });

    runComparison<Dumux::DiffMethod::numeric>("TPFA numeric experimental assembler comparison");
    runComparison<Dumux::DiffMethod::analytic>("TPFA analytic experimental assembler comparison");

    return 0;
}
