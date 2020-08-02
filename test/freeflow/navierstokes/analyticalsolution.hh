// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NavierStokesTests
 * \copydoc Dumux::NavierStokesAnalyticalSolution
 */
#ifndef DUMUX_TEST_ANALYTICAL_SOLUTION_HH
#define DUMUX_TEST_ANALYTICAL_SOLUTION_HH

#include <array>
#include <vector>
#include <type_traits>

namespace Dumux {

//! Returns a vector of the analytical velocity solution at the center of each element.
template<class Problem>
auto getVelocityAnalyticalSolution(const Problem& problem)
{
    static_assert(Problem::isMomentumProblem());
    using GridGeometry = std::decay_t<decltype(problem.gridGeometry())>;
    using GlobalPosition = typename GridGeometry::LocalView::SubControlVolumeFace::GlobalPosition;
    using PrimaryVariables = std::decay_t<decltype(problem.analyticalSolution(GlobalPosition(0.0)))>;
    using Scalar = typename PrimaryVariables::value_type;

    std::vector<std::vector<Scalar>> result(problem.gridGeometry().gridView().size(0), std::vector<Scalar>(PrimaryVariables::size()));

    for (const auto& element : elements(problem.gridGeometry().gridView()))
    {
        const auto eIdx = problem.gridGeometry().elementMapper().index(element);
        const auto center = element.geometry().center();

        // the vtk writer expects std::vector so we can't pass the PrimaryVariables (fixed-size container) directly.
        const auto sol = problem.analyticalSolution(center);
        for (int i = 0; i < PrimaryVariables::size(); ++i)
            result[eIdx][i] = sol[i];
    }

    return result;
}

//! Returns an array of vectors containing the analytical solution for a scalar value at the center of each element.
template<class Problem>
auto getScalarAnalyticalSolution(const Problem& problem)
{
    static_assert(!Problem::isMomentumProblem());
    using GridGeometry = std::decay_t<decltype(problem.gridGeometry())>;
    using GlobalPosition = typename GridGeometry::LocalView::SubControlVolumeFace::GlobalPosition;
    using PrimaryVariables = std::decay_t<decltype(problem.analyticalSolution(GlobalPosition(0.0)))>;
    using Scalar = typename PrimaryVariables::value_type;

    std::array<std::vector<Scalar>, PrimaryVariables::size()> result;
    for (auto& r : result)
        r.resize(problem.gridGeometry().gridView().size(0));

    for (const auto& element : elements(problem.gridGeometry().gridView()))
    {
        auto fvGeometry = localView(problem.gridGeometry());
        fvGeometry.bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto sol = problem.analyticalSolution(scv.dofPosition());
            for (int i = 0; i < PrimaryVariables::size(); ++i)
                result[i][scv.dofIndex()] = sol[i];
        }
    }

    return result;
}

} // end namespace Dumux

#endif
