// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test that the multi-stage assembler assembles the momentum balance of the
 *        face-centered staggered scheme like the instationary assembler
 */
#include <config.h>

#include <iostream>
#include <memory>
#include <cmath>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <test/assembly/weightedassembly.hh>

#include "properties_momentum.hh"

#if PERIODIC_DOF_MAP
#include <test/assembly/periodicdofmap.hh>

namespace Dumux::Properties {

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = Dumux::Test::GridGeometryWithPeriodicDofMap<typename GridGeometry<TypeTag, TTag::DISCRETIZATION_MODEL>::type>; };

} // end namespace Dumux::Properties
#endif

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG_MOMENTUM;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using GridManager = Dumux::GridManager<GetPropType<TypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
#if PERIODIC_DOF_MAP
    const std::size_t numPeriodicRows = Test::addInteriorPeriodicDofs(*gridGeometry);
#else
    const std::size_t numPeriodicRows = 0;
#endif
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs()), xOld(gridGeometry->numDofs());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        xOld[i] = 0.3 + 0.2*std::sin(0.61*i);
        x[i] = xOld[i] + 0.1*std::cos(0.37*i);
    }

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(xOld);

    const auto result = Test::compareWeightedAssembly<TypeTag>(problem, gridGeometry, gridVariables, x, xOld, 0.05);

    std::cout << "compared rows: " << result.numCompared << "\n"
              << "residual difference: " << result.residual << "\n"
              << "Jacobian difference: " << result.jacobian << "\n"
              << "constraint rows: " << result.numConstraintRows << " (periodic " << result.numPeriodicRows << "), difference " << result.constraint << "\n"
              << "residual without and with Jacobian: " << result.residualOnly << std::endl;

    if (result.numCompared == 0 || result.numConstraintRows == 0 || result.numPeriodicRows != numPeriodicRows || result.constraint > 1e-12 || result.residual > 1e-10 || result.jacobian > 1e-5 || result.residualOnly > 1e-12)
    {
        std::cerr << "The weighted assembly differs from the instationary assembly" << std::endl;
        return 1;
    }

    return 0;
}
