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
 * \brief Helper functions to compute L2-norm
 */
#ifndef DUMUX_TEST_L2_NORM_HH
#define DUMUX_TEST_L2_NORM_HH

#include <cmath>
#include <type_traits>
#include <dune/geometry/quadraturerules.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux {

template<class Scalar>
struct L2Norm
{
    //! Calculate the L2-Norm of the solution and hmax of the grid
    template<class Problem, class Solution>
    static Scalar computeErrorNorm(const Problem &problem, const Solution& sol, int order)
    {
        // calculate the L2-Norm and hmax
        Scalar norm = 0.0;

        // iterate over all elements
        const auto& gg = problem.fvGridGeometry();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto geometry = element.geometry();
            const auto center = geometry.center();
            const auto elemSol = elementSolution(element, sol, gg);

            // maybe exclude some elements from the norm
            static const bool excludeInnerBulk = getParam<bool>("Problem.NormExcludeInnerBulk");
            static const Scalar radius = getParam<Scalar>("SpatialParams.Radius") + 0.01;
            using GridView = std::decay_t<decltype(gg.gridView())>;
            if (int(GridView::dimension) == 3 && excludeInnerBulk && std::sqrt(center[0]*center[0] + center[1]*center[1]) < radius)
                continue;

            const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
            for(auto&& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const Scalar pe = problem.exactSolution(globalPos);
                const Scalar p = evalSolution(element, geometry, gg, elemSol, globalPos)[0];
                norm += (p - pe)*(p - pe)*qp.weight()*geometry.integrationElement(qp.position());
            }
        }

        return std::sqrt(norm);
    }

    template<class Problem>
    static Scalar computeNormalization(const Problem &problem, int order)
    {
        // calculate the L2-Norm of the exact solution vector
        Scalar norm = 0.0;

        // iterate over all elements
        const auto& gg = problem.fvGridGeometry();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto geometry = element.geometry();
            const auto center = geometry.center();

            // maybe exclude some elements from the norm
            static const bool excludeInnerBulk = getParam<bool>("Problem.NormExcludeInnerBulk");
            static const Scalar radius = getParam<Scalar>("SpatialParams.Radius") + 0.01;
            using GridView = std::decay_t<decltype(gg.gridView())>;
            if (int(GridView::dimension) == 3 && excludeInnerBulk && std::sqrt(center[0]*center[0] + center[1]*center[1]) < radius)
                continue;

            const auto& quad = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(geometry.type(), order);
            for(auto&& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const Scalar pe = problem.exactSolution(globalPos);
                norm += pe*pe*qp.weight()*geometry.integrationElement(qp.position());
            }
        }

        return std::sqrt(norm);
    }
};

} //end namespace Dumux

#endif
