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
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 */
#ifndef DUMUX_INCOMPRESSIBLE_STOKES_SOLVER_HH
#define DUMUX_INCOMPRESSIBLE_STOKES_SOLVER_HH

#include <type_traits>
#include <memory>
#include <tuple>

#include <dune/common/parametertree.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <dumux/common/math.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/symmetrizedirichlet.hh>

#include "tpfalaplaceproblem.hh"

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 * \note Uses IncompressibleStokesPreconditioner as preconditioner
 */
template<class Matrix, class Vector, class PressureGG, class LinearSolverTraits>
class IncompressibleStokesSolver
: public LinearSolver
{
    using Preconditioner = IncompressibleStokesPreconditioner<Matrix, Vector, Vector>;
public:
    /*!
     * \brief Constructor
     * \param pGridGeometry grid geometry of the pressure discretization
     * \param params a parameter tree passed along to the iterative solvers
     *
     * The parameter tree requires the (constant) viscosity in the key "Component.LiquidDynamicViscosity"
     * and the (constant) density in the key "Component.LiquidDensity"
     */
    IncompressibleStokesSolver(std::shared_ptr<const PressureGG> pGridGeometry, const Vector& dirichletDofs)
    : LinearSolver(), pGridGeometry_(std::move(pGridGeometry)), dirichletDofs_(dirichletDofs)
    {
        params_ = LinearSolverParameters<LinearSolverTraits>::createParameterTree(this->paramGroup());
        density_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDensity");
        viscosity_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDynamicViscosity");
        weight_ = getParamFromGroup<double>(this->paramGroup(), "LinearSolver.Preconditioner.MassMatrixWeight", 1.0);
        solverType_ = getParamFromGroup<std::string>(this->paramGroup(), "LinearSolver.Type", "gmres");
    }

    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        auto bTmp = b;
        auto ATmp = A;

        // make Dirichlet boundary conditions symmetric
        if (getParamFromGroup<bool>(this->paramGroup(), "LinearSolver.SymmetrizeDirichlet", false))
            symmetrizeDirichlet(ATmp, bTmp, dirichletDofs_);

        using namespace Dune::Indices;
        ATmp[_1] *= -1.0/density_;
        bTmp[_1] *= -1.0/density_;

        if (getParamFromGroup<bool>(this->paramGroup(), "LinearSolver.CheckSymmetry", false))
            std::cout << "Matrix is symmetric? "
                      << std::boolalpha << matrixIsSymmetric(ATmp) << std::endl;

        return applyIterativeSolver_(ATmp, x, bTmp);
    }

    std::string name() const
    {
        return "Block-preconditioned incompressible Stokes solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    bool applyIterativeSolver_(const Matrix& A, Vector& x, const Vector& b)
    {
        auto op = std::make_shared<Dumux::ParallelMatrixAdapter<Matrix, Vector, Vector>>(A);
        auto pop = makeTpfaLaplaceOperator<typename Preconditioner::PressureLinearOperator>(pGridGeometry_->gridView());
        auto pop2 = makePressureLinearOperator_<typename Preconditioner::PressureLinearOperator>();

        auto preconditioner = std::make_shared<Preconditioner>(op, pop, pop2, params_.sub("preconditioner"));

        std::unique_ptr<Dune::InverseOperator<Vector, Vector>> solver;
        if (solverType_ == "minres")
            solver = std::make_unique<Dune::MINRESSolver<Vector>>(op, preconditioner, params_);
        else if (solverType_ == "bicgstab")
            solver = std::make_unique<Dune::BiCGSTABSolver<Vector>>(op, preconditioner, params_);
        else if (solverType_ == "cg")
            solver = std::make_unique<Dune::CGSolver<Vector>>(op, preconditioner, params_);
        else
            solver = std::make_unique<Dune::RestartedGMResSolver<Vector>>(op, preconditioner, params_);

        auto bTmp = b;
        solver->apply(x, bTmp, result_);

        return result_.converged;
    }

    template<class LinearOperator>
    std::shared_ptr<LinearOperator> makePressureLinearOperator_()
    {
        using M = typename LinearOperator::matrix_type;

        auto massMatrix = std::make_shared<M>();
        massMatrix->setBuildMode(M::random);
        const auto numDofs = pGridGeometry_->numDofs();

        Dune::MatrixIndexSet pattern;
        pattern.resize(numDofs, numDofs);
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
            pattern.add(globalI, globalI);
        pattern.exportIdx(*massMatrix);

        const auto& gv = pGridGeometry_->gridView();
        for (const auto& element : elements(gv))
        {
            const auto eIdx = pGridGeometry_->elementMapper().index(element);
            (*massMatrix)[eIdx][eIdx] = weight_*element.geometry().volume()/viscosity_;
        }

        return std::make_shared<LinearOperator>(massMatrix);
    }

    double density_, viscosity_, weight_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::shared_ptr<const PressureGG> pGridGeometry_;
    const Vector& dirichletDofs_;
    std::string solverType_;
};

} // end namespace Dumux

#endif
