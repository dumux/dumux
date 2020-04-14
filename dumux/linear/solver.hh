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
 * \brief Base class for linear solvers
 */
#ifndef DUMUX_LINEAR_SOLVER_HH
#define DUMUX_LINEAR_SOLVER_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Base class for linear solvers
 */
class LinearSolver
{
public:
    //! export Scalar type (might be needed to set parameters from output)
    //! TODO: Do we need this?
    using Scalar = double;

    /*!
     * \brief Contruct the solver
     * \note Read parameters from the parameter tree
     *       - LinearSolver.Verbosity the verbosity level of the linear solver
     *       - LinearSolver.MaxIterations the maximum iterations of the solver
     *       - LinearSolver.ResidualReduction the residual reduction threshold, i.e. stopping criterion
     *       - LinearSolver.Preconditioner.Relaxation precondition relaxation
     *       - LinearSolver.Preconditioner.Iterations the number of preconditioner iterations
     *       - LinearSolver.Preconditioner.Verbosity the preconditioner verbosity level
     */
    LinearSolver(const std::string& paramGroup = "")
    : paramGroup_(paramGroup)
    {
        verbosity_ = getParamFromGroup<int>(paramGroup, "LinearSolver.Verbosity", 0);
        maxIter_ = getParamFromGroup<int>(paramGroup, "LinearSolver.MaxIterations", 250);
        residReduction_ = getParamFromGroup<double>(paramGroup, "LinearSolver.ResidualReduction", 1e-13);
        relaxation_ = getParamFromGroup<double>(paramGroup, "LinearSolver.Preconditioner.Relaxation", 1);
        precondIter_ = getParamFromGroup<int>(paramGroup, "LinearSolver.Preconditioner.Iterations", 1);
        precondVerbosity_ = getParamFromGroup<int>(paramGroup, "LinearSolver.Preconditioner.Verbosity", 0);
    }

    /*!
     * \brief Solve the linear system Ax = b
     * \note This has to be overloaded by the actual solver
     */
    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        DUNE_THROW(Dune::NotImplemented, "Linear solver doesn't implement a solve method!");
    }

    //! the name of the linear solver
    std::string name() const
    { return "unknown solver"; }

    //! the parameter group for getting parameter from the parameter tree
    const std::string& paramGroup() const
    { return paramGroup_; }

    //! the verbosity level
    int verbosity() const
    { return verbosity_; }

    //! set the verbosity level
    void setVerbosity(int v)
    { verbosity_ = v; }

    //! the maximum number of linear solver iterations
    int maxIter() const
    { return maxIter_; }

    //! set the maximum number of linear solver iterations
    void setMaxIter(int i)
    { maxIter_ = i; }

    //! the linear solver residual reduction
    double residReduction() const
    { return residReduction_; }

    //! set the linear solver residual reduction
    void setResidualReduction(double r)
    { residReduction_ = r; }

    //! the linear solver relaxation factor
    double relaxation() const
    { return relaxation_; }

    //! set the linear solver relaxation factor
    void setRelaxation(double r)
    { relaxation_ = r; }

    //! the number of preconditioner iterations
    int precondIter() const
    { return precondIter_; }

    //! set the number of preconditioner iterations
    void setPrecondIter(int i)
    { precondIter_ = i; }

    //! the preconditioner verbosity
    int precondVerbosity() const
    { return precondVerbosity_; }

    //! set the preconditioner verbosity
    void setPrecondVerbosity(int verbosityLevel)
    { precondVerbosity_ = verbosityLevel; }

private:
    int verbosity_;
    int maxIter_;
    double residReduction_;
    double relaxation_;
    int precondIter_;
    int precondVerbosity_;
    const std::string paramGroup_;
};

} // end namespace Dumux

#endif
