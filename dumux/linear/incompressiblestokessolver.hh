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

#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <dumux/linear/solver.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Preconditioner for the incompressible Stokes problem
 *
 *      | A   B^T |         | A    0   |           | A    0   |           | A    B^T  |
 * M =  |         |    P1 = |          |   or P2 = |          |   or P3 = |           |
 *      | B   0   |         | 0  µ^-1I |           | B  µ^-1I |           | B  µ^-1I  |
 *
 * M is the system matrix and inv(P1) is a symmetric block-diagonal preconditioner (for MINRes)
 * inv(P2) a block triangular preconditioner (for GMRes) (enabled with preconditioner.UseBlockTriangular = true)
 * and inv(P3) a symmetric block triangular preconditioner (for MINRes) (enabled with preconditioner.UseFullSchurComplement = true)
 *
 * A is inverted using an AMG (default),
 * or with a direct solver if preconditioner.UseDirectVelocitySolver is true.
 *
 * Both preconditioners guarantee bounded conditions numbers if A is inverted exactly,
 * see e.g. https://doi.org/10.1002/nla.716
 */
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

    IncompressibleStokesPreconditioner(const Matrix& A,
                                       const Vector& massMatrix,
                                       const Vector& dirichletConstraints,
                                       const Dune::ParameterTree& params,
                                       double density, double viscosity)
    : A_(A), massMatrix_(massMatrix), params_(params), density_(density), viscosity_(viscosity)
    {
        using namespace Dune::Indices;

        useDirectVelocitySolver_ = params.get<bool>("UseDirectVelocitySolver", false);
        useBlockTriangularPreconditioner_ = params.get<bool>("UseBlockTriangular", false);
        useFullSchurComplementPreconditioner_ = params.get<bool>("UseFullSchurComplement", false);

        if (useDirectVelocitySolver_)
        {
#if HAVE_UMFPACK
            using DSVel = Dune::UMFPack<std::decay_t<decltype(A[_0][_0])>>;
            const auto velVerbosity = params.get<int>("Verbosity", 0);
            directVelOp_ = std::make_unique<DSVel>(A[_0][_0], velVerbosity);
            using Precond = Dune::InverseOperator2Preconditioner<Dune::InverseOperator<SubVector<0>, SubVector<0>>>;
            velPre_ = std::make_unique<Precond>(*directVelOp_);
#else
            DUNE_THROW(Dune::Exception, "Direct solver requested but UMFPackBackend has not been found on your system.");
#endif
        }
        else
        {
            auto op = std::make_shared<LinearOperator>(A_[_0][_0]);
            velPre_ = std::make_unique<AMGSolverVelocity>(op, params);
        }
    }

    void pre (Vector&, Vector&) override {}

    void apply (Vector& v, const Vector& d) override
    {
        if (useFullSchurComplementPreconditioner_)
        {
            v = applyRightBlock_(d);
            v = applyDiagBlock_(v);
            v = applyLeftBlock_(v);
        }
        else if (useBlockTriangularPreconditioner_)
        {
            v = applyRightBlock_(d);
            v = applyDiagBlock_(v);
        }
        else
            v = applyDiagBlock_(d);
    }

    void post (Vector&) override {}

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:

    Vector applyRightBlock_(const Vector& d)
    {
        using namespace Dune::Indices;

        // right bit of Schur complement preconditioner
        // applied to d
        //
        // | I         0 |
        // | -C*inv(A) I |
        //

        auto dTmp0 = d[_0];
        auto vTmp = d; vTmp = 0.0;
        // invert velocity block (apply to first block of d)
        velPre_->pre(vTmp[_0], dTmp0);
        velPre_->apply(vTmp[_0], dTmp0);
        velPre_->post(vTmp[_0]);
        // then multiply with C
        A_[_1][_0].mv(vTmp[_0], vTmp[_1]);

        // and subtract from d
        auto v = d;
        v[_0] = vTmp[_0]; // already do A^-1 d of the diagonal block because we already computed it here
        v[_1] -= vTmp[_1];
        return v;
    }

    Vector applyDiagBlock_(const Vector& d)
    {
        using namespace Dune::Indices;

        // diagonal bit of Schur complement decomposition
        auto dTmp = d;
        auto v = d; v = 0.0;

        // invert velocity block
        if (!useFullSchurComplementPreconditioner_ && !useBlockTriangularPreconditioner_)
        {
            velPre_->pre(v[_0], dTmp[_0]);
            velPre_->apply(v[_0], dTmp[_0]);
            velPre_->post(v[_0]);
        }
        // reuse the already computed A^-1 d (see applyRightBlock_)
        else
        {
            v[_0] = dTmp[_0];
        }

        // invert pressure block (1/mu*I)^{-1}
        v[_1] = dTmp[_1];
        for (std::size_t i = 0; i < v[_1].size(); ++i)
            v[_1][i] *= viscosity_/massMatrix_[_1][i][0];

        return v;
    }

    Vector applyLeftBlock_(const Vector& d)
    {
        using namespace Dune::Indices;

        // left bit of Schur complement preconditioner
        // applied to d
        //
        // | I  -inv(A)*B |
        // | 0       I    |
        //

        auto dTmp = d;
        auto vTmp = d; vTmp = 0.0;
        // multiply with B
        A_[_0][_1].umv(dTmp[_1], vTmp[_0]);

        // invert velocity block (apply to first block of d)
        auto vTmp0 = vTmp[_0]; vTmp0 = 0.0;
        velPre_->pre(vTmp0, vTmp[_0]);
        velPre_->apply(vTmp0, vTmp[_0]);
        velPre_->post(vTmp0);

        // and subtract from d
        auto v = d;
        v[_0] -= vTmp0;
        return v;
    }

    Matrix A_;
    const Vector& massMatrix_;
    std::unique_ptr<Dune::Preconditioner<SubVector<0>, SubVector<0>>> velPre_;
    std::unique_ptr<Dune::InverseOperator<SubVector<0>, SubVector<0>>> directVelOp_;
    bool useDirectVelocitySolver_;
    bool useFullSchurComplementPreconditioner_, useBlockTriangularPreconditioner_;
    const Dune::ParameterTree& params_;
    double density_, viscosity_;
};

/*!
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 * \note Uses IncompressibleStokesPreconditioner as preconditioner
 * \note Currently only sequential
 */
template<class Matrix, class Vector>
class IncompressibleStokesSolver
: public LinearSolver
{
public:
    /*!
     * \brief Constructor
     * \param massMatrix the mass matrix as vector (for finite volume schemes this is the cell volume of the dof)
     * \param dirichletDofs vector that marks the Dirichlet dofs with a 1.0 and all other dofs with 0.0
     * \param params a parameter tree passed along to the iterative solvers
     *
     * The parameter tree requires the (constant) viscosity in the key "Component.LiquidKinematicViscosity"
     * and the (constant) density in the key "Component.LiquidDensity"
     */
    IncompressibleStokesSolver(const Vector& massMatrix, const Vector& dirichletDofs, const Dune::ParameterTree& params)
    : massMatrix_(massMatrix), dirichletDofs_(dirichletDofs), params_(params)
    {
        density_ = params.get<double>("Component.LiquidDensity");
        viscosity_ = params.get<double>("Component.LiquidKinematicViscosity")*density_;
        useBlockTriangularPreconditioner_ = params.get<bool>("preconditioner.UseBlockTriangular", false);
        useFullSchurComplementPreconditioner_ = params.get<bool>("preconditioner.UseFullSchurComplement", false);
    }

    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        auto ATmp = A;
        auto bTmp = b;

        using namespace Dune::Indices;
        using namespace Dune::Hybrid;

        // make Dirichlet boundary conditions symmetric
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
                        // value can be eliminated by a Dirichlet row
                        if (col.index() != row.index() && dirichletDofs_[j][col.index()][0] > 0.5)
                        {
                            bTmp[i][row.index()] -= bTmp[j][col.index()]*M[row.index()][col.index()][0][0];
                            M[row.index()][col.index()][0][0] = 0.0;
                        }
                    }
                }
            });
        });

        // make symmetric on big scale
        ATmp[_1][_0] *= -1.0/density_;
        ATmp[_1][_1] *= -1.0/density_;
        bTmp[_1] *= -1.0/density_;

        return applyIterativeSolver_(ATmp, x, bTmp);
    }

    std::string name() const
    {
        return "Block-preconditioned Incompressible Stokes Solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    bool applyIterativeSolver_(const Matrix& A, Vector& x, const Vector& b)
    {
        auto preconditioner = std::make_shared<IncompressibleStokesPreconditioner<Matrix, Vector>>(
            A, massMatrix_, dirichletDofs_, params_.sub("preconditioner"), density_, viscosity_
        );

        auto op = std::make_shared<Dune::MatrixAdapter<Matrix, Vector, Vector>>(A);
        std::unique_ptr<Dune::InverseOperator<Vector, Vector>> solver;

        if (useBlockTriangularPreconditioner_ && !useFullSchurComplementPreconditioner_)
            solver = std::make_unique<Dune::RestartedGMResSolver<Vector>>(op, preconditioner, params_);
        else
            solver = std::make_unique<Dune::MINRESSolver<Vector>>(op, preconditioner, params_);

        auto bTmp = b;
        Dune::InverseOperatorResult result;
        solver->apply(x, bTmp, result);

        return result.converged;
    }

    bool useFullSchurComplementPreconditioner_, useBlockTriangularPreconditioner_;
    double density_, viscosity_;
    Dune::InverseOperatorResult result_;
    const Vector massMatrix_;
    const Vector dirichletDofs_;
    const Dune::ParameterTree params_;
};

/*!
 * \ingroup Linear
 * \brief Helper function to compute the input for the IncompressibleStokesSolver
 * \param massMatrix resized vector
 * \param dirichletDofs resized vector
 * \param gridView the grid view on which we solve
 */
template<class V, class GV, class MassGG, class MomentumGG, class MomentumProblem, std::size_t i, std::size_t j>
std::tuple<V, V> computeIncompressibleStokesMassMatrixAndDirichletDofs(
    const GV& gridView,
    const MassGG& massGridGeometry, Dune::index_constant<i> freeFlowMassIndex,
    const MomentumGG& momentumGridGeometry, const MomentumProblem& momentumProblem, Dune::index_constant<j> freeFlowMomentumIndex
)
{
    // mark Dirichlet dofs and compute mass matrix
    V massMatrix;
    massMatrix[freeFlowMassIndex].resize(massGridGeometry.numDofs());
    massMatrix[freeFlowMomentumIndex].resize(momentumGridGeometry.numDofs());
    V dirichletDofs;
    dirichletDofs[freeFlowMassIndex].resize(massGridGeometry.numDofs());
    dirichletDofs[freeFlowMomentumIndex].resize(momentumGridGeometry.numDofs());
    massMatrix = 0.0;
    dirichletDofs = 0.0;

    for (const auto& element : elements(gridView))
    {
        const auto eIdx = massGridGeometry.elementMapper().index(element);
        massMatrix[freeFlowMassIndex][eIdx] = element.geometry().volume();

        auto fvGeometry = localView(momentumGridGeometry);
        fvGeometry.bind(element);
        for (const auto& scv : scvs(fvGeometry))
            massMatrix[freeFlowMomentumIndex][scv.dofIndex()] += scv.volume();

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.isFrontal() && scvf.boundary())
            {
                const auto bcTypes = momentumProblem.boundaryTypes(element, scvf);
                if (bcTypes.hasDirichlet())
                    dirichletDofs[freeFlowMomentumIndex][fvGeometry.scv(scvf.insideScvIdx()).dofIndex()] = 1.0;
            }
        }
    }

    return std::make_tuple(massMatrix, dirichletDofs);
}

} // end namespace Dumux

#endif
