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
 * \brief Iterative solver for the incompressible Stokes equations
 */
#ifndef DUMUX_LINEAR_INCOMPRESSIBLE_STOKES_SOLVER_HH
#define DUMUX_LINEAR_INCOMPRESSIBLE_STOKES_SOLVER_HH

#include <type_traits>
#include <tuple>
#include <utility>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/io.hh>

#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/timer.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/utility.hh>

#include <dumux/linear/solver.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/matrixconverter.hh>

namespace Dumux {

template<class Matrix, class Vector>
class IncompressibleStokesPreconditioner : public Dune::Preconditioner<Vector, Vector>
{
    template<std::size_t i>
    using SubVector = std::decay_t<decltype(std::declval<Vector>()[Dune::index_constant<i>{}])>;

    using A = std::decay_t<decltype(std::declval<Matrix>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<Vector>()[Dune::Indices::_0])>;

    using Comm = Dune::Amg::SequentialInformation;
    using LinearOperator = Dune::MatrixAdapter<A, U, U>;
    using Smoother = Dune::SeqSSOR<A, U, U>;
    using AMGSolverVelocity = Dune::Amg::AMG<LinearOperator, U, Smoother, Comm>;

public:
    using domain_type = Vector;
    using range_type = Vector;
    using field_type = typename Vector::field_type;
    using scalar_field_type = Dune::Simd::Scalar<field_type>;

    IncompressibleStokesPreconditioner(const Matrix& A, const Vector& massMatrix, const Vector& dirichletConstraints, const Dune::ParameterTree& params)
    : A_(A), massMatrix_(massMatrix)
    {
        using namespace Dune::Indices;

        const auto density = getParam<double>("Component.LiquidDensity");
        viscosity_ = getParam<double>("Component.LiquidKinematicViscosity")*density;

        const bool directVelSolver = getParam<bool>("LinearSolver.Velocity.UseDirectSolver", false);
        if (directVelSolver)
        {
            using DSVel = Dune::UMFPack<std::decay_t<decltype(A[_0][_0])>>;
            const auto velVerbosity = getParam<int>("LinearSolver.Velocity.Verbosity", 0);
            directVelOp_ = std::make_unique<DSVel>(A[_0][_0], velVerbosity);
            using Precond = Dune::InverseOperator2Preconditioner<Dune::InverseOperator<SubVector<0>, SubVector<0>>>;
            velOp_ = std::make_unique<Precond>(*directVelOp_);
        }
        else
        {
            // if we have pure Neumann BC add an identity operator to make the solution x to A[0][0]x = b unique
            if (dirichletConstraints.two_norm2() < 1e-7)
                for (auto rowIt = A_[_0][_0].begin(), rowEndIt = A_[_0][_0].end(); rowIt != rowEndIt; ++rowIt)
                    for (int i = 0; i < U::block_type::size(); ++i)
                        A_[_0][_0][rowIt.index()][rowIt.index()][i][i] += viscosity_*massMatrix_[_0][rowIt.index()][i];

            auto op = std::make_shared<LinearOperator>(A_[_0][_0]);
            velOp_ = std::make_unique<AMGSolverVelocity>(op, params);
        }
    }

    void pre (Vector&, Vector&) override {}

    void apply (Vector& v, const Vector& d) override
    {
        // invert velocity block
        using namespace Dune::Indices;

        auto dTmp = d[_0];
        velOp_->pre(v[_0], dTmp);
        velOp_->apply(v[_0], dTmp);
        velOp_->post(v[_0]);

        v[_1] = d[_1];
        for (std::size_t i = 0; i < v[_1].size(); ++i)
            v[_1][i] *= viscosity_/massMatrix_[_1][i][0];
    }

    void post (Vector&) override {}

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    Matrix A_;
    const Vector& massMatrix_;
    std::unique_ptr<Dune::Preconditioner<SubVector<0>, SubVector<0>>> velOp_;
    std::unique_ptr<Dune::InverseOperator<SubVector<0>, SubVector<0>>> directVelOp_;
    double viscosity_;
};

template<class Matrix, class Vector>
class IncompressibleStokesSolver
: public LinearSolver
{
    template<std::size_t i>
    using SubVector = std::decay_t<decltype(std::declval<Vector>()[Dune::index_constant<i>{}])>;

    template<std::size_t i>
    static constexpr auto blockSize = SubVector<i>::block_type::size();
public:
    template<class CouplingManager>
    IncompressibleStokesSolver(const CouplingManager& couplingManager)
    : LinearSolver()
    {
        density_ = getParam<double>("Component.LiquidDensity");
        symmetrizeDirichlet_ = getParam<bool>("LinearSolver.SymmetrizeDirichlet", true);
        useDirectSolver_ = getParam<bool>("LinearSolver.UseDirectSolver", false);

        constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
        constexpr auto massIdx = CouplingManager::freeFlowMassIndex;

        // compute Dirichlet dofs and mass matrix
        const auto& massGridGeometry = couplingManager.problem(massIdx).gridGeometry();
        const auto& momentumProblem = couplingManager.problem(momentumIdx);
        const auto& momentumGridGeometry = momentumProblem.gridGeometry();

        massMatrix_[massIdx].resize(massGridGeometry.numDofs());
        massMatrix_[momentumIdx].resize(momentumGridGeometry.numDofs());
        massMatrix_ = 0.0;
        dirichletDofs_ = massMatrix_;

        auto fvGeometryMass = localView(massGridGeometry);
        auto fvGeometryMomentum = localView(momentumGridGeometry);
        for (const auto& element : elements(massGridGeometry.gridView()))
        {
            // TODO this ignores the extrusion factor
            fvGeometryMass.bind(element);
            for (const auto& scv : scvs(fvGeometryMass))
                massMatrix_[massIdx][scv.dofIndex()] += scv.volume();

            fvGeometryMomentum.bind(element);
            for (const auto& scv : scvs(fvGeometryMomentum))
            {
                massMatrix_[momentumIdx][scv.dofIndex()] += scv.volume();

                // mark Dirichlet boundaries
                if (scv.boundary())
                {
                    auto bcTypes = momentumProblem.boundaryTypesAtPos(scv.dofPosition());
                    if (bcTypes.hasDirichlet())
                        dirichletDofs_[momentumIdx][scv.dofIndex()] = 1.0;
                }
            }
        }

        // load parameters
        using MassGridGeometry = std::decay_t<decltype(massGridGeometry)>;
        params_ = LinearSolverParameters<LinearSolverTraits<MassGridGeometry>>::createParameterTree();
    }

    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        auto ATmp = A;
        auto bTmp = b;

        using namespace Dune::Hybrid;
        using namespace Dune::Indices;

        if (bTmp[_1].size() < 17)
        {
            std::cout << "cells: " << bTmp[_1].size() << std::endl;
            std::cout << "before symmetrization:" << std::endl;
            printOutMatrix_(ATmp);
            printOutResidual_(bTmp);
        }

        // make Dirichlet boundary conditions symmetric
        if (symmetrizeDirichlet_)
        {
            forEach(std::make_index_sequence<Matrix::M()>{}, [&](const auto i)
            {
                forEach(std::make_index_sequence<Matrix::N()>{}, [&](const auto j)
                {
                    auto& M = ATmp[i][j];
                    const auto rowEnd = M.end();
                    for (auto row = M.begin(); row != rowEnd; ++row)
                    {
                        const auto colEnd = row->end();
                        for (auto col = row->begin(); col != colEnd; ++col)
                        {
                            for (int m = 0; m < blockSize<i>; ++m)
                            {
                                for (int k = 0; k < blockSize<j>; ++k)
                                {
                                    // value can be eliminated by a Dirichlet row
                                    if ((i != j || col.index() != row.index() || m != k) && dirichletDofs_[j][col.index()][k] > 0.5)
                                    {
                                        bTmp[i][row.index()][m] -= bTmp[j][col.index()][k]*M[row.index()][col.index()][m][k];
                                        M[row.index()][col.index()][m][k] = 0.0;
                                    }
                                }
                            }
                        }
                    }
                });
            });
        }

        // make symmetric on big scale
        ATmp[_1][_0] *= -1.0/density_;
        ATmp[_1][_1] *= -1.0/density_;
        bTmp[_1] *= -1.0/density_;

        // check if the matrix is symmetric
        const bool checkSymmetry = getParam<bool>("LinearSolver.CheckSymmetry", false);
        const double zeroThreshold = getParam<double>("LinearSolver.CheckSymmetryZeroThreshold", 0.0);
        if (checkSymmetry)
        {
            bool isSymmetric = true;

            using namespace Dune::Hybrid;
            forEach(std::make_index_sequence<Matrix::N()>{}, [&](const auto i)
            {
                const auto& diagBlockA = ATmp[i][i];
                for (auto rowIt = diagBlockA.begin(), rowEndIt = diagBlockA.end(); rowIt != rowEndIt; ++rowIt)
                    for (auto colIt = rowIt->begin(), colEndIt = rowIt->end(); colIt != colEndIt; ++colIt)
                        if (diagBlockA.exists(colIt.index(), rowIt.index()))
                            for (int m = 0; m < blockSize<i>; ++m)
                                for (int k = 0; k < blockSize<i>; ++k)
                                    if (Dune::FloatCmp::ne(diagBlockA[rowIt.index()][colIt.index()][m][k], diagBlockA[colIt.index()][rowIt.index()][k][m], 1e-7))
                                        if (diagBlockA[rowIt.index()][colIt.index()][m][k] > zeroThreshold && diagBlockA[colIt.index()][rowIt.index()][k][m] > zeroThreshold)
                                        {
                                            std::cout << Fmt::format("diagBlockA[{}][{}][{}][{}][{}][{}]: ", i, i, rowIt.index(), colIt.index(), m, k)
                                                      << diagBlockA[rowIt.index()][colIt.index()][m][k] << " <--> "
                                                      << diagBlockA[colIt.index()][rowIt.index()][k][m]
                                                      << std::endl;
                                            isSymmetric = false;
                                        }

                forEach(std::make_index_sequence<Matrix::N()>{}, [&](const auto j)
                {
                    if constexpr (i != j && j > i)
                    {
                        const auto& offDiag = ATmp[i][j];
                        const auto& offDiagT = ATmp[j][i];
                        for (auto rowIt = offDiag.begin(), rowEndIt = offDiag.end(); rowIt != rowEndIt; ++rowIt)
                            for (auto colIt = rowIt->begin(), colEndIt = rowIt->end(); colIt != colEndIt; ++colIt)
                                if (offDiagT.exists(colIt.index(), rowIt.index()))
                                    for (int m = 0; m < blockSize<i>; ++m)
                                        for (int k = 0; k < blockSize<j>; ++k)
                                            if (Dune::FloatCmp::ne(offDiag[rowIt.index()][colIt.index()][m][k], offDiagT[colIt.index()][rowIt.index()][k][m], 1e-7))
                                                if (offDiag[rowIt.index()][colIt.index()][m][k] > zeroThreshold && offDiagT[colIt.index()][rowIt.index()][k][m] > zeroThreshold)
                                                {
                                                    std::cout << Fmt::format("offDiag[{}][{}][{}][{}][{}][{}]: ", i, j, rowIt.index(), colIt.index(), m, k)
                                                              << offDiag[rowIt.index()][colIt.index()][m][k] << " <--> "
                                                              << offDiagT[colIt.index()][rowIt.index()][k][m]
                                                              << std::endl;
                                                    isSymmetric = false;
                                                }
                    }
                });
            });

            std::cout << Fmt::format("Matrix is symmetric? --> {}\n", isSymmetric);

            if (bTmp[_1].size() < 17)
            {
                std::cout << "cells: " << bTmp[_1].size() << std::endl;
                std::cout << "after symmetrization:" << std::endl;
                printOutMatrix_(ATmp);
                printOutResidual_(bTmp);
            }
        }

        if (useDirectSolver_)
            return applyDirectSolver_(ATmp, x, bTmp);
        else
            return applyIterativeSolver_(ATmp, x, bTmp);
    }

    std::string name() const
    {
        return "block-preconditioned Stokes solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    bool applyIterativeSolver_(const Matrix& A, Vector& x, const Vector& b)
    {
        auto op = std::make_shared<Dune::MatrixAdapter<Matrix, Vector, Vector>>(A);
        auto precond = std::make_shared<IncompressibleStokesPreconditioner<Matrix, Vector>>(
            A, massMatrix_, dirichletDofs_, params_.sub("preconditioner")
        );

        const bool symmetricA = getParam<bool>("LinearSolver.SymmetricMatrix", true);
        std::shared_ptr<Dune::InverseOperator<Vector, Vector>> solver;
        if (symmetricA)
            solver = std::make_shared<Dune::MINRESSolver<Vector>>(op, precond, params_);
        else
            solver = std::make_shared<Dune::RestartedGMResSolver<Vector>>(op, precond, params_);

        auto bTmp = b;
        solver->apply(x, bTmp, result_);

        return result_.converged;
    }

    bool applyDirectSolver_(const Matrix& A, Vector& x, const Vector& b)
    {
        Dune::Timer timer;
        const auto M = MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A);
        Dune::UMFPack<std::decay_t<decltype(M)>> solver(M, this->verbosity() > 0);

        auto bTmp = VectorConverter<Vector>::multiTypeToBlockVector(b);
        auto xTmp = VectorConverter<Vector>::multiTypeToBlockVector(x);
        std::cout << "-- Solver setup time: " << timer.elapsed() << std::endl;
        timer.reset();
        solver.apply(xTmp, bTmp, result_);

        if (result_.converged)
            VectorConverter<Vector>::retrieveValues(x, xTmp);
        std::cout << "-- Solver solve time: " << timer.elapsed() << std::endl;

        return result_.converged;
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

    double density_;
    Dune::InverseOperatorResult result_;
    Vector massMatrix_;
    Vector dirichletDofs_;
    Dune::ParameterTree params_;
    bool symmetrizeDirichlet_;
    bool useDirectSolver_;
};

} // end namespace Dumux

#endif
