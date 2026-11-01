// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test that the multi-stage multi-domain assembler assembles the coupled Navier-Stokes
 *        mass and momentum balances like the instationary multi-domain assembler, including
 *        the coupling blocks and the constraint rows
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <tuple>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/multistagemultidomainassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>

#include <test/assembly/weightedassembly_multidomain.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using MassTypeTag = Properties::TTag::TYPETAG_MASS;
    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto momentumGridGeometry = std::make_shared<GetPropType<MomentumTypeTag, Properties::GridGeometry>>(leafGridView);
    auto massGridGeometry = std::make_shared<GetPropType<MassTypeTag, Properties::GridGeometry>>(leafGridView);

    auto couplingManager = std::make_shared<CouplingManager>();
    auto momentumProblem = std::make_shared<GetPropType<MomentumTypeTag, Properties::Problem>>(momentumGridGeometry, couplingManager);
    auto massProblem = std::make_shared<GetPropType<MassTypeTag, Properties::Problem>>(massGridGeometry, couplingManager);

    const Scalar dt = 0.05;
    momentumProblem->updateTime(dt);
    massProblem->updateTime(dt);

    SolutionVector x, xOld;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    xOld = x;
    for (std::size_t i = 0; i < x[momentumIdx].size(); ++i)
    {
        xOld[momentumIdx][i] = 0.3 + 0.2*std::sin(0.61*i);
        x[momentumIdx][i] = xOld[momentumIdx][i] + 0.1*std::cos(0.37*i);
    }
    for (std::size_t i = 0; i < x[massIdx].size(); ++i)
    {
        xOld[massIdx][i] = 1.0 + 0.2*std::sin(0.83*i);
        x[massIdx][i] = xOld[massIdx][i] + 0.1*std::cos(0.29*i);
    }

    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    massGridVariables->init(xOld[massIdx]);
    momentumGridVariables->init(xOld[momentumIdx]);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 10.0*dt, false);
    using OldAssembler = Experimental::MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    OldAssembler oldAssembler(
        std::make_tuple(momentumProblem, massProblem),
        std::make_tuple(momentumGridGeometry, massGridGeometry),
        std::make_tuple(momentumGridVariables, massGridVariables),
        couplingManager, timeLoop, xOld
    );

    using MultiStageAssembler = Experimental::MultiStageMultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    MultiStageAssembler msAssembler(
        std::make_tuple(momentumProblem, massProblem),
        std::make_tuple(momentumGridGeometry, massGridGeometry),
        std::make_tuple(momentumGridVariables, massGridVariables),
        couplingManager, std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>(), xOld
    );

    const auto result = Test::compareWeightedMultiDomainAssembly(oldAssembler, msAssembler, x, xOld, dt);

    std::cout << "compared rows: " << result.numCompared << "\n"
              << "residual difference: " << result.residual << "\n"
              << "Jacobian difference: " << result.jacobian << "\n"
              << "coupling block difference: " << result.coupling << " (" << result.numCouplingEntries << " nonzero entries)\n"
              << "constraint rows: " << result.numConstraintRows << ", difference " << result.constraint << std::endl;

    if (result.numCompared == 0 || result.numConstraintRows == 0 || result.numCouplingEntries == 0
        || result.constraint > 1e-12 || result.residual > 1e-10 || result.jacobian > 1e-6)
    {
        std::cerr << "The weighted assembly differs from the instationary assembly" << std::endl;
        return 1;
    }

    return 0;
}
