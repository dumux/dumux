// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ERRORS_CVFE_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_ERRORS_CVFE_HH

#include <cmath>
#include <tuple>

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

template<class Problem, class GridVariables, class SolutionVector>
std::tuple<double, Dune::FieldVector<double, 2>> calculateL2AndH1Errors(const Problem& problem,
                                                                        const GridVariables& gridVariables,
                                                                        const SolutionVector& x,
                                                                        int order = 5)
{
    using GridGeometry = typename GridVariables::GridGeometry;
    using Extrusion = Extrusion_t<GridGeometry>;
    double totalVolume = 0.0;
    Dune::FieldVector<double, 2> errors(0.0);
    const auto& gg = problem.gridGeometry();
    auto fvGeometry = localView(gg);
    auto elemVolVars = localView(gridVariables.curGridVolVars());
    for (const auto& element : elements(gg.gridView()))
    {
        fvGeometry.bind(element);
        const auto geometry = fvGeometry.elementGeometry();

        elemVolVars.bind(element, fvGeometry, x);
        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const auto& quad = Dune::QuadratureRules<double, GridGeometry::GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            const auto& localPos = qp.position();
            const auto& qpVolumeWeight = qp.weight() * Extrusion::integrationElement(geometry, localPos);
            totalVolume += qpVolumeWeight;

            const auto& globalPos = geometry.global(localPos);
            const auto analyticalSolution = problem.analyticalSolution(globalPos);
            const auto numericalSolution = evalSolutionAtLocalPos(element, geometry, gg, elemSol, localPos);
            const auto solDiff = numericalSolution - analyticalSolution;
            errors[0] += (solDiff * solDiff) * qpVolumeWeight;

            const auto gradAnalyticalSolution = problem.gradAnalyticalSolution(globalPos);
            const auto gradNumericalSolution = evalGradientsAtLocalPos(element, geometry, gg, elemSol, localPos);
            const auto gradDiff = gradNumericalSolution - gradAnalyticalSolution;
            double gradDiffNorm = 0.0;
            for (int i = 0; i < gradDiff.size(); ++i)
                gradDiffNorm += gradDiff[i] * gradDiff[i];

            errors[1] += (solDiff * solDiff + gradDiffNorm) * qpVolumeWeight;
        }
    }
    errors[0] = std::sqrt(errors[0]);
    errors[1] = std::sqrt(errors[1]);

    return {totalVolume, errors};
}

} // end namespace Dumux

#endif
