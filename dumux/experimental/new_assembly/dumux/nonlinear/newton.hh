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
 * \ingroup NonLinear
 * \copydoc Dumux::NewtonSolver
 */
#ifndef DUMUX_NONLINEAR_NEWTON_HH
#define DUMUX_NONLINEAR_NEWTON_HH

#include <cmath>
#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/variables.hh>
#include <dumux/experimental/new_assembly/dumux/common/pdesolver.hh>

#include <dumux/experimental/new_assembly/dumux/linear/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

#ifndef DOXYGEN
// forward declaration
namespace Dune {

template<typename... Rows>
class MultiTypeBlockVector;

} // namespace Dune
#endif // DOXYGEN


namespace Dumux {

// TODO: Remove this header and use <dumux/nonlinear/newton.hh> instead!
// Just a temporary implementation to check the minimum requirements.

#ifndef DOXYGEN
namespace Detail {

// actual newton has this specialized/generalized
// should be generalized s.t. it works for any supported vector type
auto maxRelativeShift(const double& v1, const double& v2)
{
    using std::abs;
    using std::max;
    return abs(v1 - v2)/max(1.0, abs(v1 + v2)*0.5);
}

template<typename V1, typename V2>
auto maxRelativeShift(const V1& v1, const V2& v2)
{ return maxRelativeShift(v1[0], v2[0]); }

template<typename... Rows>
auto maxRelativeShift(const Dune::MultiTypeBlockVector<Rows...>& v1,
                      const Dune::MultiTypeBlockVector<Rows...>& v2)
{
    using V = Dune::MultiTypeBlockVector<Rows...>;
    using Scalar = typename V::field_type;

    using std::max;
    Scalar shift = 0.0;
    Dune::Hybrid::forEach(std::make_index_sequence<V::N()>{}, [&] (auto i) {
        shift = max(shift, maxRelativeShift(v1[i], v2[i]));
    });

    return shift;
}

} // namespace Detail
#endif // DOXYGEN


namespace Concepts {

template<typename Solver, typename F>
concept NewtonCompatibleLinearSolver
    = LinearizableFunction<F>
    and Variables<typename F::Domain>
    and LinearSolver<Solver,
                     Variables::DofsType<typename F::Domain>,
                     typename F::Linearization::Value,
                     typename F::Linearization::Derivative>;

} // namespace Concepts


/*!
 * \ingroup NonLinear
 * \brief Implementation of a Newton solver
 */
template<Concepts::LinearizableFunction F,
         Concepts::NewtonCompatibleLinearSolver<F> LinearSolver>
class NewtonSolver : public PDESolver<typename F::Domain>
{
public:
    using Variables = typename F::Domain;
    using Scalar = Dumux::Variables::ScalarType<Variables>;

    struct Options
    {
        Scalar maxRelativeShift = 1e-8;
        int maxIterations = 20;
    };

    NewtonSolver(std::shared_ptr<F> function,
                 std::shared_ptr<LinearSolver> solver,
                 Options opts = Options{})
    : function_(std::move(function))
    , solver_(std::move(solver))
    , opts_(std::move(opts))
    {}

private:
    bool solve_(Variables& vars) override
    {
        bool converged = false;
        std::size_t iterationCount = 0;

        auto uLastIter = Dumux::Variables::getDofs(vars);
        auto deltaU = Dumux::Variables::getDofs(vars);
        auto uNew = Dumux::Variables::getDofs(vars);

        while (!converged && iterationCount < opts_.maxIterations)
        {
            LinearSystem::fill(deltaU, 0.0);
            if (!solveLinearSystem_(function_->linearizeAt(vars), deltaU))
            {
                std::cout << "Linear solver did not converge" << std::endl;
                return false;
            }

            LinearSystem::subtract(uNew, deltaU);
            Dumux::Variables::update(vars, uNew);

            iterationCount++;
            const auto maxShift = Detail::maxRelativeShift(uNew, uLastIter);
            converged = maxShift < opts_.maxRelativeShift;
            std::cout << "Newton iteration " << iterationCount << " done. "
                      << "Max relative shift = " << maxShift
                      << std::endl;
            uLastIter = uNew;
        }

        return converged;
    }

    template<Concepts::Linearization L, typename DofVector>
    bool solveLinearSystem_(const L& linearization, DofVector& deltaU) const
    {
        solver_->setLinearOperator(linearization.derivative());
        return solver_->solve(deltaU, linearization.value());
    }

    std::shared_ptr<F> function_;
    std::shared_ptr<LinearSolver> solver_;
    Options opts_;
};

} // namespace Dumux

#endif
