// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \brief Checks that the RASPEN Jacobian is the derivative of the RASPEN residual.
 *
 * A central finite difference of \f$ F_{RASPEN} \f$ along a fixed direction is compared against
 * \f$ J_{RASPEN} \f$ applied to the same direction.
 *
 * The agreement cannot be driven to zero by refining the step: under DiffMethod::numeric the assembled
 * \f$ J \f$ is itself a finite-difference approximation carrying a relative error of order \f$ 10^{-7} \f$,
 * so the comparison plateaus at the accuracy of the reference rather than converging at second order.
 * What the plateau level does establish is whether the operator is the right one, so the test is stated
 * as agreement against two controls that must be markedly worse: omitting the fringe coupling entirely
 * (leaving the identity, which is what a fringe-blind implementation computes), and loosening the local
 * solves so that the Jacobian is cached at an unconverged iterate.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/nonlinear/preconditioning/partition.hh>
#include <dumux/nonlinear/preconditioning/restrictedadditiveschwarzoperator.hh>

#include "2p/properties.hh"

namespace {

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << what << std::endl;
        ++failures;
    }
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TwoPIncompressibleTpfa;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector u(gridGeometry->numDofs());
    problem->applyInitialSolution(u);

    // a genuinely nonlinear but physical state: saturations stay well inside the mobile range, since
    // driving them past the residual saturation makes the relative permeabilities, and with them the
    // local Jacobian blocks, degenerate
    for (std::size_t i = 0; i < u.size(); ++i)
    {
        u[i][0] += 1e3*std::sin(3.0*double(i));
        u[i][1] = 0.3 + 0.2*std::sin(7.0*double(i) + 1.0);
    }

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(u);

    using Index = NonlinearPreconditioning::Partition::Index;
    const auto numElements = gridGeometry->gridView().size(0);

    std::vector<std::vector<Index>> influence(numElements);
    for (Index i = 0; i < static_cast<Index>(numElements); ++i)
        for (const auto& dataJ : gridGeometry->connectivityMap()[i])
            influence[i].push_back(dataJ.globalJ);

    std::vector<NonlinearPreconditioning::SubdomainIndex> cores(numElements);
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        const auto eIdx = gridGeometry->elementMapper().index(element);
        cores[eIdx] = element.geometry().center()[0] < 3.0 ? 0 : 1;
    }

    auto partition = std::make_shared<NonlinearPreconditioning::Partition>(influence, cores, 1);

    // the subdomain problems must be transient: without a storage term the incompressible two-phase
    // system does not determine the saturation pointwise and the local blocks are singular
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, 250.0, 3000.0);
    SolutionVector uOld(u);

    NonlinearPreconditioning::RestrictedAdditiveSchwarzOperator<TypeTag, DiffMethod::numeric> rasOperator(
        problem, gridGeometry, gridVariables, partition, influence);
    rasOperator.setTimeLoop(timeLoop);
    rasOperator.setPreviousSolution(uOld);
    rasOperator.setLocalTolerance(1e-10);
    rasOperator.setLocalMaxIterations(50);

    // a direction with both pressure and saturation content
    SolutionVector direction(u);
    for (std::size_t i = 0; i < direction.size(); ++i)
    {
        direction[i][0] = 1e3*std::cos(2.0*double(i) + 0.3);
        direction[i][1] = 1e-2*std::cos(5.0*double(i) + 0.7);
    }

    SolutionVector base(u);
    check(rasOperator.evaluate(u, base), "all subdomain solves converged at the base point");

    SolutionVector exact(u);
    rasOperator.applyJacobian(direction, exact);

    const double h = 5e-3;
    SolutionVector reference(u);
    {
        SolutionVector plus(u), minus(u);
        for (std::size_t i = 0; i < u.size(); ++i)
        {
            plus[i].axpy(h, direction[i]);
            minus[i].axpy(-h, direction[i]);
        }

        SolutionVector fPlus(u), fMinus(u);
        const bool okPlus = rasOperator.evaluate(plus, fPlus);
        const bool okMinus = rasOperator.evaluate(minus, fMinus);
        check(okPlus && okMinus, "perturbed subdomain solves converged");

        for (std::size_t i = 0; i < u.size(); ++i)
            for (int eq = 0; eq < 2; ++eq)
                reference[i][eq] = (fPlus[i][eq] - fMinus[i][eq])/(2.0*h);
    }

    const auto relativeError = [&] (const SolutionVector& candidate)
    {
        double error = 0.0, scale = 0.0;
        for (std::size_t i = 0; i < u.size(); ++i)
            for (int eq = 0; eq < 2; ++eq)
            {
                error = std::max(error, std::abs(candidate[i][eq] - reference[i][eq]));
                scale = std::max(scale, std::abs(reference[i][eq]));
            }
        return error/scale;
    };

    const auto errorExact = relativeError(exact);

    // the fringe coupling is the entire content of the RASPEN Jacobian beyond the identity
    SolutionVector noFringe(direction);
    const auto errorNoFringe = relativeError(noFringe);

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  relative error, RASPEN Jacobian      = " << errorExact << std::endl;
    std::cout << "  relative error, fringe term omitted  = " << errorNoFringe << std::endl;

    check(errorExact < 1e-5, "RASPEN Jacobian agrees with the finite difference");
    check(errorNoFringe > 100.0*errorExact,
          "omitting the fringe coupling is markedly worse, so the test discriminates");

    // caching the Jacobian at an unconverged local iterate must also degrade the agreement
    rasOperator.setLocalTolerance(1e-1);
    rasOperator.setLocalMaxIterations(1);
    SolutionVector loose(u), looseJacobian(u);
    rasOperator.evaluate(u, loose);
    rasOperator.applyJacobian(direction, looseJacobian);
    const auto errorLoose = relativeError(looseJacobian);
    std::cout << "  relative error, unconverged local solves = " << errorLoose << std::endl;
    check(errorLoose > 10.0*errorExact,
          "an unconverged local solve degrades the Jacobian, so the base point matters");

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "RASPEN Jacobian matches the finite-difference derivative" << std::endl;
    return 0;
}
