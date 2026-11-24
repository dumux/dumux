// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test that the experimental assembler also supports FC staggered models.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/fvassembler.hh>

#include "../../freeflow/navierstokes/donea/properties_momentum.hh"

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

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init([](Dune::ParameterTree& params)
    {
        params["Problem.Name"] = "experimental_assembler_fcstaggered";
        params["Problem.EnableGravity"] = "false";
        params["Problem.EnableInertiaTerms"] = "false";
        params["Problem.UseNeumann"] = "false";
        params["Component.LiquidDensity"] = "1.0";
        params["Component.LiquidKinematicViscosity"] = "1.0";
        params["Assembly.NumericDifference.BaseEpsilon"] = "1e-1";
        params["SpatialParams.Temperature"] = "298.0";
    });

    using TypeTag = Dumux::Properties::TTag::DoneaTestMomentum;
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
    for (std::size_t i = 0; i < x.size(); ++i)
        x[i][0] = 0.1 + 0.03*std::sin(1.7*static_cast<double>(i))
                      - 0.02*std::cos(0.9*static_cast<double>(i));

    auto referenceGridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    auto experimentalGridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    referenceGridVariables->init(x);
    experimentalGridVariables->init(x);

    using ReferenceAssembler = Dumux::FVAssembler<TypeTag, Dumux::DiffMethod::numeric>;
    using ExperimentalAssembler = Dumux::Experimental::Assembler<TypeTag, Dumux::DiffMethod::numeric>;

    ReferenceAssembler reference(problem, gridGeometry, referenceGridVariables);
    ExperimentalAssembler experimental(problem, gridGeometry, experimentalGridVariables);
    reference.setLinearSystem();
    experimental.setLinearSystem();

    reference.assembleJacobianAndResidual(x);
    experimental.assembleJacobianAndResidual(x);

    const auto jacobianError = maxAbsDiff(reference.jacobian(), experimental.jacobian());
    const auto residualError = maxAbsDiffVector(reference.residual(), experimental.residual());

    std::cout << "FC staggered numeric experimental assembler comparison: jacobian "
              << jacobianError << ", residual " << residualError << std::endl;

    if (jacobianError > 1e-20 || residualError > 1e-20)
        DUNE_THROW(Dune::InvalidStateException, "Experimental::Assembler differs from FVAssembler for FCStaggered");

    return 0;
}
