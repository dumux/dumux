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
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include "tpfalaplaceproblem.hh"

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 * \note Uses IncompressibleStokesPreconditioner as preconditioner
 */
template<class Matrix, class Vector, class VelocityGG, class PressureGG>
class IncompressibleStokesSolver
: public LinearSolver
{
    // "[Scalar and matrix-vector products] are easy to compute in a non-overlapping distributed setting,
    // if matrix and residuals are kept in _additive [unique]_ representation,
    // and iterates and directions are kept in _consistent_ representation." (Sander 2020 Dune)
    // "The consistent representation of a vector is obtained from its additive representation by taking the latter,
    // and for each vertex[/face] in the border partition, add the corresponding entries from the other subdomains" (Sander 2020 Dune)
    using Communication = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using Preconditioner = IncompressibleStokesPreconditioner<Matrix, Vector, Vector, Communication>;
public:
    /*!
     * \brief Constructor
     * \param vGridGeometry grid geometry of the velocity discretization
     * \param pGridGeometry grid geometry of the pressure discretization
     * \param params a parameter tree passed along to the iterative solvers
     *
     * The parameter tree requires the (constant) viscosity in the key "Component.LiquidDynamicViscosity"
     * and the (constant) density in the key "Component.LiquidDensity"
     */
    IncompressibleStokesSolver(std::shared_ptr<const VelocityGG> vGridGeometry,
                               std::shared_ptr<const PressureGG> pGridGeometry,
                               const Vector& dirichletDofs)
    : LinearSolver()
    , vGridGeometry_(std::move(vGridGeometry))
    , pGridGeometry_(std::move(pGridGeometry))
    , dirichletDofs_(dirichletDofs)
    {
        params_ = LinearSolverParameters<LinearSolverTraits<VelocityGG>>::createParameterTree(this->paramGroup());
        density_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDensity");
        viscosity_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDynamicViscosity");
        weight_ = getParamFromGroup<double>(this->paramGroup(), "LinearSolver.Preconditioner.MassMatrixWeight", 1.0);
        solverType_ = getParamFromGroup<std::string>(this->paramGroup(), "LinearSolver.Type", "gmres");

        pHelperVelocity_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits<VelocityGG>>>(
            vGridGeometry_->gridView(), vGridGeometry_->dofMapper()
        );
        pHelperPressure_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits<PressureGG>>>(
            pGridGeometry_->gridView(), pGridGeometry_->dofMapper()
        );

        commVelocity_ = std::make_shared<Communication>(pHelperVelocity_->gridView().comm(), Dune::SolverCategory::nonoverlapping);
        pHelperVelocity_->createParallelIndexSet(*commVelocity_);
        commPressure_ = std::make_shared<Communication>(pHelperPressure_->gridView().comm(), Dune::SolverCategory::overlapping);
        pHelperPressure_->createParallelIndexSet(*commPressure_);

        const auto comms = std::array<std::shared_ptr<const Communication>, 2>({ commVelocity_, commPressure_ });
        scalarProduct_ = std::make_shared<Dumux::ParallelScalarProduct<Vector, Communication, Vector::size()>>(comms);
    }

    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        auto bTmp = b;
        auto ATmp = A;

        return applyIterativeSolver_(ATmp, x, bTmp);
    }

    Scalar norm(const Vector& b) const
    {
        return scalarProduct_->norm(b);
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
    bool applyIterativeSolver_(Matrix& A, Vector& x, Vector& b)
    {
        // make Dirichlet boundary conditions symmetric
        // do this first because its easier if the operator is still purely local
        // otherwise we might sum up wrong contributions at border face dofs adjacent to the domain boundary
        if (getParamFromGroup<bool>(this->paramGroup(), "LinearSolver.SymmetrizeDirichlet", false))
            symmetrizeDirichlet(A, b, dirichletDofs_);

        // make matrix and right-hand side consistent
        {
            using namespace Dune::Indices;
            using M = std::decay_t<decltype(std::declval<Matrix>()[_0][_0])>;
            using U = std::decay_t<decltype(std::declval<Vector>()[_0])>;
            using PTraits = typename LinearSolverTraits<VelocityGG>::template ParallelNonoverlapping<M, U>;
            prepareLinearAlgebraParallel<LinearSolverTraits<VelocityGG>, PTraits>(A[_0][_0], b[_0], *pHelperVelocity_);
        }

        // make Matrix symmetric
        using namespace Dune::Indices;
        A[_1] *= -1.0/density_;
        b[_1] *= -1.0/density_;

        const auto comms = std::array<std::shared_ptr<const Communication>, 2>({ commVelocity_, commPressure_ });
        auto op = std::make_shared<Dumux::ParallelStokesOperator<Matrix, Vector, Vector, Communication>>(A, comms);

        auto pop = makeTpfaLaplaceOperator<typename Preconditioner::PressureLinearOperator>(pGridGeometry_->gridView(), *commPressure_);
        auto pop2 = makePressureLinearOperator_<typename Preconditioner::PressureLinearOperator>(*commPressure_);
        auto preconditioner = std::make_shared<Preconditioner>(op, pop, pop2, comms, params_.sub("preconditioner"));

        params_["verbose"] = pGridGeometry_->gridView().comm().rank() == 0 ? params_["verbose"] : "0";
        std::unique_ptr<Dune::InverseOperator<Vector, Vector>> solver;
        if (solverType_ == "minres")
            solver = std::make_unique<Dune::MINRESSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else if (solverType_ == "bicgstab")
            solver = std::make_unique<Dune::BiCGSTABSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else if (solverType_ == "cg")
            solver = std::make_unique<Dune::CGSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else
            solver = std::make_unique<Dune::RestartedGMResSolver<Vector>>(op, scalarProduct_, preconditioner, params_);

        solver->apply(x, b, result_);

        return result_.converged;
    }

    template<class LinearOperator, class C>
    std::shared_ptr<LinearOperator> makePressureLinearOperator_(const C& comm)
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
            // do not modify ghosts
            if (element.partitionType() == Dune::GhostEntity)
                (*massMatrix)[eIdx][eIdx] = 1.0;
            else
                (*massMatrix)[eIdx][eIdx] = weight_*element.geometry().volume()/viscosity_;
        }

        return std::make_shared<LinearOperator>(massMatrix, comm);
    }

    void printOutMatrix_(const Matrix& A) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Matrix::N()>(), [&](const auto i) {
            forEach(std::make_index_sequence<Matrix::M()>(), [&](const auto j) {
                if (A[i][j].nonzeroes() > 0)
                    Dune::printmatrix(std::cout, A[i][j], "A" + std::to_string(i()) + std::to_string(j()), "");
                else
                    std::cout << "A" << i() << j() << ": 0" << std::endl;
            });
        });
    }

    void printOutResidual_(const Vector& b) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Vector::size()>(), [&](const auto i) {
            Dune::printvector(std::cout, b[i], "b" + std::to_string(i()), "");
        });
    }

    double density_, viscosity_, weight_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::shared_ptr<const VelocityGG> vGridGeometry_;
    std::shared_ptr<const PressureGG> pGridGeometry_;
    const Vector& dirichletDofs_;
    std::string solverType_;

#if HAVE_MPI
    std::unique_ptr<ParallelISTLHelper<LinearSolverTraits<VelocityGG>>> pHelperVelocity_;
    std::unique_ptr<ParallelISTLHelper<LinearSolverTraits<PressureGG>>> pHelperPressure_;

    std::shared_ptr<Communication> commVelocity_;
    std::shared_ptr<Communication> commPressure_;
#endif
    std::shared_ptr<Dune::ScalarProduct<Vector>> scalarProduct_;
};

} // end namespace Dumux

#endif
