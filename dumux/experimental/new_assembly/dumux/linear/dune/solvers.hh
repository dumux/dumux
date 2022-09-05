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
 * \ingroup Linear
 * \brief Linear solvers based on Dune istl.
 */
#ifndef DUMUX_LINEAR_DUNE_SOLVERS_HH
#define DUMUX_LINEAR_DUNE_SOLVERS_HH

#include <cmath>
#include <utility>

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioner.hh>

// include traits specializations for dune linear systems
#include <dumux/experimental/new_assembly/dumux/linear/dune/systemtraits.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/operatortraits.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

namespace Dumux {

//! \todo: implement solver using istl solver factory
template<typename Domain, typename Range>
class MatrixFreeDuneIterativeSolver;

//! \todo: implement solver using istl solver factory
template<typename Domain, typename Range>
class DuneIterativeSolver;

//! \todo: remove once generic factory-based solver is available
template<typename D>
class DuneCGSolver
{
    using Scalar = LinearSystem::ScalarType<D>;

    class NullPreconditioner : public Dune::Preconditioner<D, D>
    {
    public:
        void pre(D& x, D& b) override {}
        void apply(D& x, const D& r) override { x = r; }
        void post(D& x) override {}
        Dune::SolverCategory::Category category() const override
        { return Dune::SolverCategory::Category::sequential; }
    };

public:
    using Domain = D;
    using Range = Domain;
    using LinearOperator = Dune::LinearOperator<Domain, Range>;

    struct Options
    {
        Scalar reduction = 1e-13;
        int maxIterations = 250;
        bool verbosity = 0;
    };

    DuneCGSolver(const LinearOperator& op,
                 Options opts = Options{})
    : opts_(std::move(opts))
    , op_(&op)
    {}

    DuneCGSolver(Options opts = Options{})
    : opts_(std::move(opts))
    {}

    bool solve(Domain& u, const Range& r) const
    {
        auto rhs = r;
        NullPreconditioner prec;
        Dune::InverseOperatorResult result;
        Dune::CGSolver<Domain> solver(
            *op_,
            prec,
            opts_.reduction,
            opts_.maxIterations,
            opts_.verbosity
        );
        solver.apply(u, rhs, result);
        return result.converged;
    }

    auto scalarProduct(const Range& a, const Range& b) const
    {
        using std::sqrt;
        return sqrt(a.dot(b));
    }

    void setLinearOperator(const LinearOperator& op)
    { op_ = &op; }

private:
    Options opts_;
    const LinearOperator* op_{nullptr};
};

} // namespace Dumux

#endif
