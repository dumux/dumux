// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Dumux sequential linear solver backends
 */
#ifndef DUMUX_SEQ_SOLVER_BACKEND_HH
#define DUMUX_SEQ_SOLVER_BACKEND_HH

#include <type_traits>
#include <tuple>
#include <utility>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelmatrixadapter.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A general solver backend allowing arbitrary preconditioners and solvers.
 *
 * This class is used as a base class for specific solver-preconditioner
 * combinations. Several parameters from the group LinearSolver are read to
 * customize the solver and preconditioner:
 *
 * - Verbosity: determines how verbose the linear solver should print output.
 * - MaxIterations: the maximum number of iterations for the linear solver.
 * - ResidualReduction: the threshold for declaration of convergence.
 * - PreconditionerRelaxation: relaxation parameter for the preconditioner.
 * - PreconditionerIterations: usually specifies the number of times the
 *                             preconditioner is applied. In case of ILU(n),
 *                             it specifies the order of the applied ILU.
 */
class IterativePreconditionedSolverImpl
{
public:

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    [[deprecated("Removed after 3.8. Use solver from istlsolvers.hh")]]
    static bool solve(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                      const std::string& modelParamGroup = "")
    {
        Preconditioner precond(A, s.precondIter(), s.relaxation());

        // make a linear operator from a matrix
        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter linearOperator(A);

        Solver solver(linearOperator, precond, s.residReduction(), s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    [[deprecated("Removed after 3.8. Use solver from istlsolvers.hh")]]
    static bool solveWithGMRes(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                               const std::string& modelParamGroup = "")
    {
        // get the restart threshold
        const int restartGMRes = getParamFromGroup<int>(modelParamGroup, "LinearSolver.GMResRestart", 10);

        Preconditioner precond(A, s.precondIter(), s.relaxation());

        // make a linear operator from a matrix
        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter linearOperator(A);

        Solver solver(linearOperator, precond, s.residReduction(), restartGMRes, s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    [[deprecated("Removed after 3.8. Use solver from istlsolvers.hh")]]
    static bool solveWithILU0Prec(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                                  const std::string& modelParamGroup = "")
    {
        Preconditioner precond(A, s.relaxation());

        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter operatorA(A);

        Solver solver(operatorA, precond, s.residReduction(), s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    // solve with RestartedGMRes (needs restartGMRes as additional argument)
    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    [[deprecated("Removed after 3.8. Use solver from istlsolvers.hh")]]
    static bool solveWithILU0PrecGMRes(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                                       const std::string& modelParamGroup = "")
    {
        // get the restart threshold
        const int restartGMRes = getParamFromGroup<int>(modelParamGroup, "LinearSolver.GMResRestart", 10);

        Preconditioner precond(A, s.relaxation());

        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter operatorA(A);

        Solver solver(operatorA, precond, s.residReduction(), restartGMRes, s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    // solve with generic parameter tree
    template<class Preconditioner, class Solver, class Matrix, class Vector>
    static bool solveWithParamTree(const Matrix& A, Vector& x, const Vector& b,
                                   const Dune::ParameterTree& params)
    {
        // make a linear operator from a matrix
        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        const auto linearOperator = std::make_shared<MatrixAdapter>(A);

        auto precond = std::make_shared<Preconditioner>(linearOperator, params.sub("preconditioner"));
        Solver solver(linearOperator, precond, params);

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }
};

/*!
 * \ingroup Linear
 * \brief Returns the block level for the preconditioner for a given matrix
 *
 * \tparam M The matrix.
 */
template<class M>
constexpr std::size_t preconditionerBlockLevel() noexcept
{
    return isMultiTypeBlockMatrix<M>::value ? 2 : 1;
}

/*!
 * \ingroup Linear
 * \brief Solver for simple block-diagonal matrices (e.g. from explicit time stepping schemes)
 *
 * Solver: Single Jacobi iteration
 * Preconditioner: Unity
 */
class ExplicitDiagonalSolver : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        Vector rhs(b);
        constexpr auto precondBlockLevel = preconditionerBlockLevel<Matrix>();
        Dune::SeqJac<Matrix, Vector, Vector, precondBlockLevel> jac(A, 1, 1.0);
        jac.pre(x, rhs);
        jac.apply(x, rhs);
        jac.post(x);
        return true;
    }

    std::string name() const
    {
        return "Explicit diagonal matrix solver";
    }
};

/*!
 * \name Solver for MultiTypeBlockMatrix's
 */
// \{

/*!
 * \ingroup Linear
 * \brief A Uzawa preconditioned BiCGSTAB solver for saddle-point problems
 */
template <class LinearSolverTraits>
class UzawaBiCGSTABBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = SeqUzawa<Matrix, Vector, Vector>;
        using Solver = Dune::BiCGSTABSolver<Vector>;
        static const auto solverParams = LinearSolverParameters<LinearSolverTraits>::createParameterTree(this->paramGroup());
        return IterativePreconditionedSolverImpl::template solveWithParamTree<Preconditioner, Solver>(A, x, b, solverParams);
    }

    std::string name() const
    {
        return "Uzawa preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioner
 */
template<class M, class X, class Y, int blockLevel = 2>
class BlockDiagILU0Preconditioner : public Dune::Preconditioner<X, Y>
{
    template<std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using BlockILU = Dune::SeqILU<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>, blockLevel-1>;

    using ILUTuple = typename makeFromIndexedType<std::tuple, BlockILU, std::make_index_sequence<M::N()> >::type;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::decay_t<M>;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param m The (multi type block) matrix to operate on
       \param w The relaxation factor
     */
    BlockDiagILU0Preconditioner(const M& m, double w = 1.0)
    : BlockDiagILU0Preconditioner(m, w, std::make_index_sequence<M::N()>{})
    {
        static_assert(blockLevel >= 2, "Only makes sense for MultiTypeBlockMatrix!");
    }

    void pre (X& v, Y& d) final {}

    void apply (X& v, const Y& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(ilu_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(ilu_).apply(v[i], d[i]);
        });
    }

    void post (X&) final {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const final
    {
        return Dune::SolverCategory::sequential;
    }

private:
    template<std::size_t... Is>
    BlockDiagILU0Preconditioner (const M& m, double w, std::index_sequence<Is...> is)
    : ilu_(std::make_tuple(BlockILU<Is>(m[Dune::index_constant<Is>{}][Dune::index_constant<Is>{}], w)...))
    {}

    ILUTuple ilu_;
};


/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioned BiCGSTABSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
class BlockDiagILU0BiCGSTABSolver : public LinearSolver
{

public:
    using LinearSolver::LinearSolver;

    template<class Matrix, class Vector>
    bool solve(const Matrix& M, Vector& x, const Vector& b)
    {
        BlockDiagILU0Preconditioner<Matrix, Vector, Vector> preconditioner(M);
        Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector> op(M);
        Dune::BiCGSTABSolver<Vector> solver(op, preconditioner, this->residReduction(),
                                            this->maxIter(), this->verbosity());
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    std::string name() const
    { return "block-diagonal ILU0-preconditioned BiCGSTAB solver"; }

private:
    Dune::InverseOperatorResult result_;
};

/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioned RestartedGMResSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
class BlockDiagILU0RestartedGMResSolver : public LinearSolver
{

public:
    using LinearSolver::LinearSolver;

    template<int precondBlockLevel = 2, class Matrix, class Vector>
    bool solve(const Matrix& M, Vector& x, const Vector& b)
    {
        BlockDiagILU0Preconditioner<Matrix, Vector, Vector> preconditioner(M);
        Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector> op(M);
        static const int restartGMRes = getParamFromGroup<int>(this->paramGroup(), "LinearSolver.GMResRestart");
        Dune::RestartedGMResSolver<Vector> solver(op, preconditioner, this->residReduction(), restartGMRes,
                                                  this->maxIter(), this->verbosity());
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    std::string name() const
    { return "block-diagonal ILU0-preconditioned restarted GMRes solver"; }

private:
    Dune::InverseOperatorResult result_;
};

/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioner
 */
template<class M, class X, class Y, int blockLevel = 2>
class BlockDiagAMGPreconditioner : public Dune::Preconditioner<X, Y>
{
    template<std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using Smoother = Dune::SeqSSOR<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template<std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template<std::size_t i>
    using ScalarProduct = Dune::SeqScalarProduct<VecBlockType<i>>;

    template<std::size_t i>
    using BlockAMG = Dune::Amg::AMG<LinearOperator<i>, VecBlockType<i>, Smoother<i>>;

    using AMGTuple = typename makeFromIndexedType<std::tuple, BlockAMG, std::make_index_sequence<M::N()> >::type;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::decay_t<M>;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param lop The linear operator
       \param c The criterion
       \param sa The smoother arguments
     */
    template<class LOP, class Criterion, class SmootherArgs>
    BlockDiagAMGPreconditioner(const LOP& lop, const Criterion& c, const SmootherArgs& sa)
    : BlockDiagAMGPreconditioner(lop, c, sa, std::make_index_sequence<M::N()>{})
    {
        static_assert(blockLevel >= 2, "Only makes sense for MultiTypeBlockMatrix!");
    }

    /*!
     * \brief Prepare the preconditioner.
     *
     * A solver solves a linear operator equation A(v)=d by applying
     * one or several steps of the preconditioner. The method pre()
     * is called before the first apply operation.
     * d and v are right hand side and solution vector of the linear
     * system respectively. It may. e.g., scale the system, allocate memory or
     * compute a (I)LU decomposition.
     * Note: The ILU decomposition could also be computed in the constructor
     * or with a separate method of the derived method if several
     * linear systems with the same matrix are to be solved.
     *
     * \param v The left hand side of the equation.
     * \param d The right hand side of the equation.
     */
    void pre (X& v, Y& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(amg_).pre(v[i], d[i]);
        });
    }

    /*!
     * \brief Apply one step of the preconditioner to the system A(v)=d.
     *
     * On entry v=0 and d=b-A(x) (although this might not be
     * computed in that way. On exit v contains the update, i.e
     * one step computes \f$ v = M^{-1} d \f$ where \f$ M \f$ is the
     * approximate inverse of the operator \f$ A \f$ characterizing
     * the preconditioner.
     * \param v The update to be computed
     * \param d The current defect.
     */
    void apply (X& v, const Y& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(amg_).apply(v[i], d[i]);
        });
    }

    /*!
     * \brief Clean up.
     *
     * This method is called after the last apply call for the
     * linear system to be solved. Memory may be deallocated safely
     * here. v is the solution of the linear equation.
     *
     * \param v The right hand side of the equation.
     */
    void post (X& v) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(amg_).post(v[i]);
        });
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const final
    {
        return Dune::SolverCategory::sequential;
    }

private:
    template<class LOP, class Criterion, class SmootherArgs, std::size_t... Is>
    BlockDiagAMGPreconditioner (const LOP& lop, const Criterion& c, const SmootherArgs& sa, std::index_sequence<Is...> is)
    : amg_(std::make_tuple(BlockAMG<Is>(*std::get<Is>(lop), *std::get<Is>(c), *std::get<Is>(sa))...))
    {}

    AMGTuple amg_;
};

/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioned BiCGSTABSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
class BlockDiagAMGBiCGSTABSolver : public LinearSolver
{
    template<class M, std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<class X, std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template<class M, class X, std::size_t i>
    using Smoother = Dune::SeqSSOR<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;

    template<class M, class X, std::size_t i>
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother<M, X, i>>::Arguments;

    template<class M, std::size_t i>
    using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<DiagBlockType<M, i>, Dune::Amg::FirstDiagonal>>;

    template<class M, class X, std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;

public:
    using LinearSolver::LinearSolver;

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<class Matrix, class Vector>
    bool solve(const Matrix& m, Vector& x, const Vector& b)
    {
        //! \todo Check whether the default accumulation mode atOnceAccu is needed.
        //! \todo make parameters changeable at runtime from input file / parameter tree
        Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(3);
        params.setDebugLevel(this->verbosity());

        auto criterion = makeCriterion_<Criterion, Matrix>(params, std::make_index_sequence<Matrix::N()>{});
        auto smootherArgs = makeSmootherArgs_<SmootherArgs, Matrix, Vector>(std::make_index_sequence<Matrix::N()>{});

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Matrix::N()>{}, [&](const auto i)
        {
            auto& args = std::get<decltype(i)::value>(smootherArgs);
            args->iterations = 1;
            args->relaxationFactor = 1;
        });

        auto linearOperator = makeLinearOperator_<LinearOperator, Matrix, Vector>(m, std::make_index_sequence<Matrix::N()>{});

        BlockDiagAMGPreconditioner<Matrix, Vector, Vector> preconditioner(linearOperator, criterion, smootherArgs);

        Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector> op(m);
        Dune::BiCGSTABSolver<Vector> solver(op, preconditioner, this->residReduction(),
                                            this->maxIter(), this->verbosity());
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    std::string name() const
    { return "block-diagonal AMG-preconditioned BiCGSTAB solver"; }

private:
    template<template<class M, std::size_t i> class Criterion, class Matrix, class Params, std::size_t... Is>
    auto makeCriterion_(const Params& p, std::index_sequence<Is...>)
    {
        return std::make_tuple(std::make_shared<Criterion<Matrix, Is>>(p)...);
    }

    template<template<class M, class X, std::size_t i> class SmootherArgs, class Matrix, class Vector, std::size_t... Is>
    auto makeSmootherArgs_(std::index_sequence<Is...>)
    {
        return std::make_tuple(std::make_shared<SmootherArgs<Matrix, Vector, Is>>()...);
    }

    template<template<class M, class X, std::size_t i> class LinearOperator, class Matrix, class Vector, std::size_t... Is>
    auto makeLinearOperator_(const Matrix& m, std::index_sequence<Is...>)
    {
        return std::make_tuple(std::make_shared<LinearOperator<Matrix, Vector, Is>>(m[Dune::index_constant<Is>{}][Dune::index_constant<Is>{}])...);
    }

    Dune::InverseOperatorResult result_;
};

// \}

} // end namespace Dumux

#endif
