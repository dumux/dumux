// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 */
#ifndef DUMUX_LINEAR_STOKES_SOLVER_HH
#define DUMUX_LINEAR_STOKES_SOLVER_HH

#include <type_traits>
#include <memory>
#include <tuple>

#include <dune/common/parametertree.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/exceptions.hh>

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <dumux/common/math.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelmatrixadapter.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/linear/symmetrize_constraints.hh>

namespace Dumux::Detail {

/*!
 * \ingroup Linear
 * \brief A Stokes preconditioner (saddle-point problem) for the problem
 * \f$
 \begin{pmatrix} A & B \\ C & 0 \end{pmatrix}
 \begin{pmatrix} u \\ p \end{pmatrix} =
 \begin{pmatrix} f \\ g \end{pmatrix},
 * \f$
 *
 * where A is the discrete velocity operator and B and C are discrete gradient and divergence operators.
 * This preconditioner is especially suited for solving the incompressible Stokes equations.
 *
 * \tparam M Type of the matrix.
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level
 */
template<class M, class X, class Y, int l = 2>
class StokesPreconditioner : public Dune::Preconditioner<X,Y>
{
    static_assert(Dumux::isMultiTypeBlockMatrix<M>::value && M::M() == 2 && M::N() == 2, "Expects a 2x2 MultiTypeBlockMatrix.");
    static_assert(l == 2, "StokesPreconditioner expects a block level of 2 for Matrix type M.");

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;

    using P = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_1][Dune::Indices::_1])>;
    using V = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>;

    enum class Mode { symmetric, triangular,  diagonal };

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief Scalar type underlying the field_type.
    using scalar_field_type = Dune::Simd::Scalar<field_type>;
    //! \brief the type of the pressure operator
    using PressureLinearOperator = Dune::MatrixAdapter<P,V,V>;

    /*!
     * \brief Constructor
     * \param fullLinearOperator the Stokes linear operator
     * \param pressureLinearOperator the linear operator to be used for the pressure space preconditioner
     * \param params a parameter tree for the preconditioner configuration
     */
    StokesPreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& fullLinearOperator,
        const std::shared_ptr<const PressureLinearOperator>& pressureLinearOperator,
        const Dune::ParameterTree& params
    )
    : matrix_(fullLinearOperator->getmat())
    , pmatrix_(pressureLinearOperator->getmat())
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    {
        const auto mode = getParamFromGroup<std::string>(
            paramGroup_, "LinearSolver.Preconditioner.Mode", "Diagonal"
        );

        if (mode == "Symmetric")
            mode_ = Mode::symmetric;
        else if (mode == "Triangular")
            mode_ = Mode::triangular;
        else
            mode_ = Mode::diagonal;

        initPreconditioner_(params);
    }

    /*!
     * \brief Prepare the preconditioner.
     */
    void pre(X& update, Y& currentDefect) override {}

    /*!
     * \brief Apply the preconditioner
     *
     * \param update The update to be computed.
     * \param currentDefect The current defect.
     *
     * The currentDefect has be be in a consistent representation,
     * Definition 2.3 Blatt and Bastian (2009) https://doi.org/10.1504/IJCSE.2008.021112
     * The update is initially zero. At exit the update has to be
     * in a consistent representation. This usually requires communication.
     */
    void apply(X& update, const Y& currentDefect) override
    {
        using namespace Dune::Indices;

        if (mode_ == Mode::symmetric)
        {
            update = applyRightBlock_(currentDefect);
            update = applyDiagBlock_(update);
            update = applyLeftBlock_(update);
        }
        else if (mode_ == Mode::triangular)
        {
            update = applyRightBlock_(currentDefect);
            update = applyDiagBlock_(update);
        }
        else
            update = applyDiagBlock_(currentDefect);
    }

    /*!
     * \brief Clean up.
     */
    void post(X& update) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    X applyRightBlock_(const Y& d)
    {
        using namespace Dune::Indices;

        // right bit of LDU decomposition
        // applied to d
        //
        // | I         0 |
        // | -C*inv(A) I |
        //

        auto dTmp0 = d[_0];
        auto vTmp = d; vTmp = 0.0;

        // invert velocity block (apply to first block of d)
        applyPreconditionerForA_(vTmp[_0], dTmp0);

        // then multiply with C
        matrix_[_1][_0].mv(vTmp[_0], vTmp[_1]);

        // and subtract from d
        auto v = d;
        v[_0] = vTmp[_0]; // already do A^-1 d of the diagonal block because we already computed it here
        v[_1] -= vTmp[_1];
        return v;
    }

    X applyDiagBlock_(const Y& d)
    {
        using namespace Dune::Indices;

        // diagonal middle bit of LDU decomposition
        auto dTmp = d;
        auto v = d; v = 0.0;

        // invert velocity block
        if (mode_ == Mode::diagonal)
            applyPreconditionerForA_(v[_0], dTmp[_0]);

        // reuse the already computed A^-1 d (see applyRightBlock_)
        else
            v[_0] = dTmp[_0];

        // invert pressure block
        applyPreconditionerForP_(v[_1], dTmp[_1]);

        return v;
    }

    X applyLeftBlock_(const Y& d)
    {
        using namespace Dune::Indices;

        // left bit of LDU decomposition
        // applied to d
        //
        // | I  -inv(A)*B |
        // | 0       I    |
        //

        auto dTmp = d;
        auto vTmp = d; vTmp = 0.0;

        // multiply with B
        matrix_[_0][_1].umv(dTmp[_1], vTmp[_0]);

        // invert velocity block (apply to first block of d)
        auto vTmp0 = vTmp[_0]; vTmp0 = 0.0;
        applyPreconditionerForA_(vTmp0, vTmp[_0]);

        // and subtract from d
        auto v = d;
        v[_0] -= vTmp0;

        return v;
    }

    void initPreconditioner_(const Dune::ParameterTree& params)
    {
        using namespace Dune::Indices;

        if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.DirectSolverForVelocity", false))
        {
#if HAVE_UMFPACK
            directSolver_ = std::make_shared<Dune::UMFPack<A>>(matrix_[_0][_0], verbosity_);
            using Wrap = Dune::InverseOperator2Preconditioner<Dune::InverseOperator<U, U>>;
            preconditionerForA_ = std::make_shared<Wrap>(*directSolver_);
#else
            DUNE_THROW(Dune::InvalidStateException, "Selected direct solver but UMFPack is not available.");
#endif
        }
        else
        {
            using VelLinearOperator = Dune::MatrixAdapter<A, U, U>;
            auto lopV = std::make_shared<VelLinearOperator>(matrix_[_0][_0]);
            preconditionerForA_ = std::make_shared<
                Dune::Amg::AMG<VelLinearOperator, U, Dumux::ParMTSSOR<A,U,U>>
            >(lopV, params);
        }

        using PressJacobi = Dumux::ParMTJac<P, V, V>;
        preconditionerForP_ = std::make_shared<PressJacobi>(pmatrix_, 1, 1.0);
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForA_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForA_->pre(sol, rhs);
        preconditionerForA_->apply(sol, rhs);
        preconditionerForA_->post(sol);
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForP_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForP_->pre(sol, rhs);
        preconditionerForP_->apply(sol, rhs);
        preconditionerForP_->post(sol);
    }

    //! \brief The matrix we operate on.
    const M& matrix_;
    //! \brief The matrix we operate on.
    const P& pmatrix_;
    //! \brief The verbosity level
    const int verbosity_;

    std::shared_ptr<Dune::Preconditioner<U, U>> preconditionerForA_;
    std::shared_ptr<Dune::Preconditioner<V, V>> preconditionerForP_;
    std::shared_ptr<Dune::InverseOperator<U, U>> directSolver_;

    const std::string paramGroup_;
    Mode mode_;
};

} // end namespace Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 * \note Uses StokesPreconditioner as preconditioner (tailored to the incompressible Stokes problem)
 * \note No MPI parallelization implemented, some shared-memory parallelism is enabled
 */
template<class Matrix, class Vector, class VelocityGG, class PressureGG>
class StokesSolver
: public LinearSolver
{
    using Preconditioner = Detail::StokesPreconditioner<Matrix, Vector, Vector>;
public:
    /*!
     * \brief Constructor
     * \param vGridGeometry grid geometry of the velocity discretization
     * \param pGridGeometry grid geometry of the pressure discretization
     * \param dirichletDofs a vector (same size and shape as right hand side) where the dirichlet dofs are marked with 1.0
     *                      and all other entries are 0.0.
     * \param paramGroup group prefix when looking up keys in the parameter tree
     */
    StokesSolver(std::shared_ptr<const VelocityGG> vGridGeometry,
                 std::shared_ptr<const PressureGG> pGridGeometry,
                 const Vector& dirichletDofs,
                 const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , vGridGeometry_(std::move(vGridGeometry))
    , pGridGeometry_(std::move(pGridGeometry))
    , dirichletDofs_(dirichletDofs)
    {
        params_ = LinearSolverParameters<LinearSolverTraits<VelocityGG>>::createParameterTree(this->paramGroup());
        density_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDensity");
        viscosity_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDynamicViscosity");
        weight_ = getParamFromGroup<double>(this->paramGroup(), "LinearSolver.Preconditioner.MassMatrixWeight", 1.0);
        solverType_ = getParamFromGroup<std::string>(this->paramGroup(), "LinearSolver.Type", "gmres");
        scalarProduct_ = std::make_shared<Dune::ScalarProduct<Vector>>();
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
        return "Block-preconditioned Stokes solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    bool applyIterativeSolver_(Matrix& A, Vector& x, Vector& b)
    {
        // make Dirichlet boundary conditions symmetric
        if (getParamFromGroup<bool>(this->paramGroup(), "LinearSolver.SymmetrizeDirichlet", true))
            symmetrizeConstraints(A, b, dirichletDofs_);

        // make Matrix symmetric on the block-scale
        using namespace Dune::Indices;
        A[_1] *= -1.0/density_;
        b[_1] *= -1.0/density_;

        auto op = std::make_shared<Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector>>(A);
        auto pop = makePressureLinearOperator_<typename Preconditioner::PressureLinearOperator>();
        auto preconditioner = std::make_shared<Preconditioner>(op, pop, params_.sub("preconditioner"));
        params_["verbose"] = pGridGeometry_->gridView().comm().rank() == 0 ? params_["verbose"] : "0";

        // defaults to restarted GMRes
        std::unique_ptr<Dune::InverseOperator<Vector, Vector>> solver;
        if (solverType_ == "minres")
            solver = std::make_unique<Dune::MINRESSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else if (solverType_ == "bicgstab")
            solver = std::make_unique<Dune::BiCGSTABSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else if (solverType_ == "gmres")
            solver = std::make_unique<Dune::RestartedGMResSolver<Vector>>(op, scalarProduct_, preconditioner, params_);
        else
            DUNE_THROW(Dune::NotImplemented, "Solver choice " << solverType_ << " is not implemented");

        solver->apply(x, b, result_);

        return result_.converged;
    }

    template<class LinearOperator>
    std::shared_ptr<LinearOperator> makePressureLinearOperator_()
    {
        using M = typename LinearOperator::matrix_type;
        auto massMatrix = createMassMatrix_<M>();
        return std::make_shared<LinearOperator>(massMatrix);
    }

    template<class M>
    std::shared_ptr<M> createMassMatrix_()
    {
        auto massMatrix = std::make_shared<M>();
        massMatrix->setBuildMode(M::random);
        const auto numDofs = pGridGeometry_->numDofs();

        Dune::MatrixIndexSet pattern;
        pattern.resize(numDofs, numDofs);
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
            pattern.add(globalI, globalI);
        pattern.exportIdx(*massMatrix);

        const auto& gv = pGridGeometry_->gridView();
        auto fvGeometry = localView(*pGridGeometry_);
        for (const auto& element : elements(gv))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                using Extrusion = Extrusion_t<PressureGG>;
                const auto dofIndex = scv.dofIndex();
                if (element.partitionType() == Dune::GhostEntity) // do not modify ghosts
                    (*massMatrix)[dofIndex][dofIndex] = 1.0;
                else
                    (*massMatrix)[dofIndex][dofIndex] += weight_*Extrusion::volume(fvGeometry, scv)/(2.0*viscosity_);
            }
        }

        return massMatrix;
    }

    double density_, viscosity_, weight_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::shared_ptr<const VelocityGG> vGridGeometry_;
    std::shared_ptr<const PressureGG> pGridGeometry_;
    const Vector& dirichletDofs_;
    std::string solverType_;

    std::shared_ptr<Dune::ScalarProduct<Vector>> scalarProduct_;
};

} // end namespace Dumux

#endif
