// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Discrete error norms of a shallow water solution with respect to an exact solution
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_ERRORS_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_ERRORS_HH

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <type_traits>

#include <dune/grid/common/partitionset.hh>

#include <dumux/io/format.hh>

namespace Dumux::ShallowWaterTest {

/*!
 * \ingroup ShallowWaterTests
 * \brief Discrete \f$ L^1 \f$ and \f$ L^2 \f$ norms of the errors in water depth and velocity
 *
 * Both norms are normalized with the area of the domain,
 * \f$ \| e \|_1 = \sum_i |V_i| \, |e_i| / \sum_i |V_i| \f$ and
 * \f$ \| e \|_2 = \left( \sum_i |V_i| \, e_i^2 / \sum_i |V_i| \right)^{1/2} \f$,
 * where \f$ e_i \f$ is the error at the center of the cell \f$ V_i \f$. The norms of the exact
 * solution, evaluated in the same way, give the relative errors.
 * The components are ordered as the primary variables \f$ (h, u, v) \f$.
 */
template<class Scalar>
struct Errors
{
    std::array<Scalar, 3> l1 = {};
    std::array<Scalar, 3> l2 = {};
    std::array<Scalar, 3> exactL1 = {};
    std::array<Scalar, 3> exactL2 = {};
    std::size_t numElements = 0;
};

/*!
 * \ingroup ShallowWaterTests
 * \brief Compute the error norms of the solution vector
 * \param gridGeometry The cell-centered finite volume grid geometry
 * \param x The solution vector with the primary variables \f$ (h, u, v) \f$
 * \param exactSolution A callable returning the exact \f$ (h, u, v) \f$ for a given element
 */
template<class GridGeometry, class SolutionVector, class ExactSolution>
auto computeErrors(const GridGeometry& gridGeometry, const SolutionVector& x, const ExactSolution& exactSolution)
{
    using Scalar = std::decay_t<decltype(x[0][0])>;
    Errors<Scalar> errors;
    Scalar area = 0.0;

    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        fvGeometry.bindElement(element);
        const auto exact = exactSolution(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            for (int i = 0; i < 3; ++i)
            {
                using std::abs;
                const auto error = abs(x[scv.dofIndex()][i] - exact[i]);
                errors.l1[i] += error*scv.volume();
                errors.l2[i] += error*error*scv.volume();
                errors.exactL1[i] += abs(exact[i])*scv.volume();
                errors.exactL2[i] += exact[i]*exact[i]*scv.volume();
            }
            area += scv.volume();
        }
        ++errors.numElements;
    }

    const auto& comm = gridGeometry.gridView().comm();
    area = comm.sum(area);
    errors.numElements = comm.sum(errors.numElements);
    for (int i = 0; i < 3; ++i)
    {
        using std::sqrt;
        errors.l1[i] = comm.sum(errors.l1[i])/area;
        errors.l2[i] = sqrt(comm.sum(errors.l2[i])/area);
        errors.exactL1[i] = comm.sum(errors.exactL1[i])/area;
        errors.exactL2[i] = sqrt(comm.sum(errors.exactL2[i])/area);
    }

    return errors;
}

/*!
 * \ingroup ShallowWaterTests
 * \brief Print the error norms of a stationary solution in the format read by the convergence test scripts
 */
template<class Scalar>
void printErrors(const Errors<Scalar>& errors)
{
    std::cout << Fmt::format("Errors (h, u, v) on {} elements:\n", errors.numElements)
              << Fmt::format("  L1: {:.8e} {:.8e} {:.8e}\n", errors.l1[0], errors.l1[1], errors.l1[2])
              << Fmt::format("  L2: {:.8e} {:.8e} {:.8e}\n", errors.l2[0], errors.l2[1], errors.l2[2])
              << "Norms of the exact solution (h, u, v):\n"
              << Fmt::format("  L1: {:.8e} {:.8e} {:.8e}\n", errors.exactL1[0], errors.exactL1[1], errors.exactL1[2])
              << Fmt::format("  L2: {:.8e} {:.8e} {:.8e}\n", errors.exactL2[0], errors.exactL2[1], errors.exactL2[2])
              << std::flush;
}

/*!
 * \ingroup ShallowWaterTests
 * \brief Print the error norms at the given time in the format read by the convergence test scripts
 */
template<class Scalar>
void printErrors(const Errors<Scalar>& errors, Scalar time)
{
    std::cout << Fmt::format("At t = {:.6g} s:\n", time);
    printErrors(errors);
}

} // end namespace Dumux::ShallowWaterTest

#endif
