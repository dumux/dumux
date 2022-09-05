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
 * \brief Linear solver implementations using Eigen3
 */
#ifndef DUMUX_LINEAR_EIGEN_SOLVERS_HH
#define DUMUX_LINEAR_EIGEN_SOLVERS_HH

#include <config.h>

#if HAVE_EIGEN3

#include <utility>
#include <cassert>

#include <Eigen/Core>
#include <Eigen/Sparse>

// include traits specializations for eigen linear systems
#include <dumux/experimental/new_assembly/dumux/linear/eigen/systemtraits.hh>
#include <dumux/experimental/new_assembly/dumux/linear/eigen/operatortraits.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

namespace Dumux {

//! CGSolver using Eigen
template<typename EigenMatrix = Eigen::SparseMatrix<double>,
         typename EigenVector = Eigen::VectorXd>
class EigenCGSolver
{
    using Solver = Eigen::ConjugateGradient<EigenMatrix>;
    using Scalar = LinearSystem::ScalarType<EigenVector>;

public:
    using Domain = EigenVector;
    using Range = EigenVector;
    using LinearOperator = const EigenMatrix;

    struct Options
    {
        Scalar tolerance = 1e-13;
        int maxIterations = 250;
    };

    EigenCGSolver(const LinearOperator& op,
                  const Options& opts = Options{})
    : op_(&op)
    , solver_(makeSolver_(opts))
    {}

    EigenCGSolver(const Options& opts = Options{})
    : op_(nullptr)
    , solver_(makeSolver_(opts))
    {}

    bool solve(Domain& u, const Range& r) const
    {
        assert(op_ && solver_);
        u = solver_->solve(r);
        return solver_->info() == Eigen::Success;
    }

    auto scalarProduct(const Range& a, const Range& b) const
    {
        using std::sqrt;
        return sqrt(a.dot(b));
    }

    void setLinearOperator(const LinearOperator& op)
    {
        op_ = &op;
        solver_->compute(*op_);
    }

private:
    auto makeSolver_(const Options& opts)
    {
        auto solver = std::make_unique<Solver>();
        solver->setTolerance(opts.tolerance);
        solver->setMaxIterations(opts.maxIterations);
        return solver;
    }

    const LinearOperator* op_{nullptr};
    std::unique_ptr<Solver> solver_{nullptr};
};

} // namespace Dumux

#endif // HAVE_EIGEN3
#endif
