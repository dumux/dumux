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
#include <dumux/linear/parallelhelpers.hh>

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


template<class Vector, class PreconditionerTuple>
class TuplePreconditioner : public Dune::Preconditioner<Vector, Vector>
{
public:
    TuplePreconditioner(const PreconditionerTuple& precTuple)
    : precTuple_(precTuple)
    {}

    void pre (Vector& v, Vector& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(v)), [&](const auto i)
        {
            std::get<i>(precTuple_)->pre(v[i], d[i]);
        });
    }

    void apply (Vector& v, const Vector& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(v)), [&](const auto i)
        {
            std::get<i>(precTuple_)->apply(v[i], d[i]);
        });
    }

    void post (Vector& v) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(v)), [&](const auto i)
        {
            std::get<i>(precTuple_)->post(v[i]);
        });
    }

    Dune::SolverCategory::Category category() const final
    {
        if (std::get<0>(precTuple_)->category() == Dune::SolverCategory::sequential)
            return Dune::SolverCategory::sequential;

        return Dune::SolverCategory::overlapping;
    }

private:
    PreconditionerTuple precTuple_;
};


template<class Vector, class Matrix, class LinearOperatorTuple>
class TupleLinearOperator: public Dune::LinearOperator<Vector, Vector> {
public:
    //! The type of the domain of the operator.
    typedef Vector domain_type;
    //! The type of the range of the operator.
    typedef Vector range_type;
    //! The field type of the operator.
    typedef typename Vector::field_type field_type;

    TupleLinearOperator (const LinearOperatorTuple& lops, const Matrix& m)
    : lops_(lops), m_(m)
    {}

    /*! \brief apply operator to x:  \f$ y = A(x) \f$
     *          The input vector is consistent and the output must also be
     *       consistent on the interior+border partition.
     */
    void apply (const Vector& x, Vector& y) const override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            std::get<i>(lops_)->apply(x[i], y[i]);

            forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto j)
            {
                if (i != j)
                    m_[i][j].umv(x[j], y[i]);
            });
        });
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd (field_type alpha, const Vector& x, Vector& y) const override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            std::get<i>(lops_)->applyscaleadd(alpha, x[i], y[i]);

            forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto j)
            {
                if (i != j)
                    m_[i][j].usmv(alpha, x[j], y[i]);
            });
        });
    }

    //! Category of the linear operator (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::overlapping;
    }

private:
    const LinearOperatorTuple& lops_;
    const Matrix& m_;
};


template<class Vector, class ScalarProductTuple>
class TupleScalarProduct : public Dune::ScalarProduct<Vector>
{
    using field_type = typename Dune::ScalarProduct<Vector>::field_type;
    using real_type = typename Dune::ScalarProduct<Vector>::real_type;

public:
    /*!
     * \param comm The communication object for syncing overlap and copy
     * data points.
     * \param cat parallel solver category (nonoverlapping or overlapping)
     */
    TupleScalarProduct (const ScalarProductTuple& sps)
    : sps_(sps)
    {}

    /*! \brief Dot product of two vectors.
     *       It is assumed that the vectors are consistent on the interior+border
     *       partition.
     */
    field_type dot (const Vector& x, const Vector& y) const override
    {
        field_type result(0);

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            result += std::get<i>(sps_)->dot(x[i], y[i]);
        });

        return result;
    }

    /*! \brief Norm of a right-hand side vector.
     *       The vector must be consistent on the interior+border partition
     */
    real_type norm (const Vector& x) const override
    {
        using std::sqrt;
        return sqrt(dot(x, x));
    }

    //! Category of the scalar product (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::overlapping;
    }

private:
    const ScalarProductTuple& sps_;
};

/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioned BiCGSTABSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
template <class LinearSolverTraitsTuple>
class BlockDiagAMGGMResSolver : public LinearSolver
{
    template<class M, std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<class X, std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using LinearSolverTraits = std::tuple_element_t<i, LinearSolverTraitsTuple>;

    template<class X, std::size_t i>
    using ScalarProduct = Dune::ScalarProduct<VecBlockType<X, i>>;

    template<std::size_t i>
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;

    template<class X, std::size_t i>
    using LinearOperator = Dune::LinearOperator<VecBlockType<X, i>, VecBlockType<X, i>>;

    template<class X, std::size_t i>
    using Preconditioner = Dune::Preconditioner<VecBlockType<X, i>, VecBlockType<X, i>>;

    template<std::size_t i>
    using ParallelHelperSP = std::shared_ptr<ParallelISTLHelper<LinearSolverTraits<i>>>;

    template<std::size_t i>
    using ParallelHelper = ParallelISTLHelper<LinearSolverTraits<i>>;

    static constexpr auto numBlocks = std::tuple_size_v<LinearSolverTraitsTuple>;
    using ParallelHelperTuple = typename makeFromIndexedType<std::tuple,
                                                             ParallelHelperSP,
                                                             std::make_index_sequence<numBlocks>
                                                            >::type;

public:
    using LinearSolver::LinearSolver;

    template<class GridView, class DofMapper, class String>
    BlockDiagAMGGMResSolver(const GridView& gridView,
                               const DofMapper& dofMapper,
                               const String& paramGroup = "")
    : BlockDiagAMGGMResSolver(gridView, dofMapper, paramGroup, std::make_index_sequence<numBlocks>{})
    {}

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<class Matrix, class Vector>
    bool solve(Matrix& m, Vector& x, Vector& b)
    {
#if HAVE_MPI
        solveSequentialOrParallel_(m, x, b);
#else
        solveSequential_(m, x, b);
#endif
        firstCall_ = false;
        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    void reset()
    {
        firstCall_ = true;
    }

    std::string name() const
    { return "block-diagonal AMG preconditioned GMRes solver"; }

private:

#if HAVE_MPI
    template <class SolverTraits, class ParallelTraits,
              class MatrixBlock, class VectorBlock, class Comm, class LOP, class SP, class PH, class Prec>
    void prepareAlgebra_(MatrixBlock& diagBlock, VectorBlock& rhsBlock, std::shared_ptr<Comm>& comm,
                         LOP& linearOperator, SP& scalarProduct,
                         PH& parHelper, Prec& preconditioner)
    {
        Dune::SolverCategory::Category category;

        if constexpr (ParallelTraits::isNonOverlapping)
        {
            using GridView = typename SolverTraits::GridView;
            using DofMapper = typename SolverTraits::DofMapper;
            static constexpr int dofCodim = SolverTraits::dofCodim;
            ParallelMatrixHelper<MatrixBlock, GridView, DofMapper, dofCodim> matrixHelper(parHelper.gridView(), parHelper.dofMapper());
            matrixHelper.extendMatrix(diagBlock, [&parHelper](auto idx){ return parHelper.isGhost(idx); });
            matrixHelper.sumEntries(diagBlock);

            ParallelVectorHelper<GridView, DofMapper, dofCodim> vectorHelper(parHelper.gridView(), parHelper.dofMapper());
            vectorHelper.makeNonOverlappingConsistent(rhsBlock);

            category = Dune::SolverCategory::nonoverlapping;
        }
        else
        {
            category = Dune::SolverCategory::overlapping;
        }

        comm = std::make_shared<Comm>(parHelper.gridView().comm(), category);
        parHelper.createParallelIndexSet(*comm);
        linearOperator = std::make_shared<typename ParallelTraits::LinearOperator>(diagBlock, *comm);
        scalarProduct = std::make_shared<typename ParallelTraits::ScalarProduct>(*comm);

        using SeqSmoother = Dune::SeqSSOR<MatrixBlock, VectorBlock, VectorBlock>;
        using Smoother = typename ParallelTraits::template Preconditioner<SeqSmoother>;
        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        SmootherArgs args;
        args.iterations = 1;
        args.relaxationFactor = 1;

        Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
        params.setDebugLevel(this->verbosity());
        params.setDefaultValuesIsotropic(SolverTraits::GridView::dimension);

        using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixBlock, Dune::Amg::FirstDiagonal>>;
        Criterion criterion(params);

        // Cast the linear operator from a pointer to the base class
        // to a pointer to the actually employed derived class.
        using ParallelLinearOperator = typename ParallelTraits::LinearOperator;
        auto lop = std::dynamic_pointer_cast<ParallelLinearOperator>(linearOperator);

        using AMG = Dune::Amg::AMG<ParallelLinearOperator, VectorBlock, Smoother, Comm>;
        preconditioner = std::make_shared<AMG>(*lop, criterion, args, *comm);
    }

    template<class Matrix, class Vector>
    void solveSequentialOrParallel_(Matrix& m, Vector& x, Vector& b)
    {
        using LSTraits = LinearSolverTraits<0>;

        if (LSTraits::canCommunicate && isParallel_)
        {
            solveParallel_(m, x, b);
        }
        else
        {
            solveSequential_(m, x, b);
        }
    }

    template<class Matrix, class Vector>
    void solveParallel_(Matrix& m, Vector& x, Vector& b)
    {
        auto linearOperators = makeLinearOperators_<LinearOperator, Vector>(std::make_index_sequence<numBlocks>{});
        auto scalarProducts = makeScalarProducts_<ScalarProduct, Vector>(std::make_index_sequence<numBlocks>{});
        auto preconditioners = makePreconditioners_<Preconditioner, Vector>(std::make_index_sequence<numBlocks>{});
        auto comms = makeComms_<Comm>(std::make_index_sequence<numBlocks>{});

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(x)), [&](const auto i)
        {
            auto& diagBlock = m[Dune::index_constant<i>{}][Dune::index_constant<i>{}];
            auto& rhsBlock = b[Dune::index_constant<i>{}];

            using LSTraits = LinearSolverTraits<i>;
            using DiagBlock = DiagBlockType<Matrix, i>;
            using RHSBlock = VecBlockType<Vector, i>;

            auto& comm = std::get<i>(comms);
            auto& linearOperator = std::get<i>(linearOperators);
            auto& scalarProduct = std::get<i>(scalarProducts);
            auto& parHelper = *std::get<i>(parHelpers_);

            if (LSTraits::isNonOverlapping(parHelper.gridView()))
            {
                using PTraits = typename LSTraits::template ParallelNonoverlapping<DiagBlock, RHSBlock>;
                prepareAlgebra_<LSTraits, PTraits>(diagBlock, rhsBlock, comm, linearOperator, scalarProduct,
                                                   parHelper, std::get<i>(preconditioners));
            }
            else
            {
                using PTraits = typename LSTraits::template ParallelOverlapping<DiagBlock, RHSBlock>;
                prepareAlgebra_<LSTraits, PTraits>(diagBlock, rhsBlock, comm, linearOperator, scalarProduct,
                                                   parHelper, std::get<i>(preconditioners));
            }
        });

        TuplePreconditioner<Vector, decltype(preconditioners)> prec(preconditioners);

        TupleLinearOperator<Vector, Matrix, decltype(linearOperators)> op(linearOperators, m);

        TupleScalarProduct<Vector, decltype(scalarProducts)> sp(scalarProducts);

        auto rank = Dune::MPIHelper::getCollectiveCommunication().rank();
        static const int restartGMRes = getParamFromGroup<double>(this->paramGroup(), "LinearSolver.GMResRestart", 10);
        Dune::RestartedGMResSolver<Vector> solver(op, sp, prec, this->residReduction(), restartGMRes,
                                                  this->maxIter(), rank == 0 ? this->verbosity() : 0);

        auto bTmp(b);
        solver.apply(x, bTmp, result_);
    }
#endif

    template<class Matrix, class Vector>
    void solveSequential_(Matrix& m, Vector& x, Vector& b)
    {
        auto linearOperators = makeLinearOperators_<LinearOperator, Vector>(std::make_index_sequence<numBlocks>{});
        auto preconditioners = makePreconditioners_<Preconditioner, Vector>(std::make_index_sequence<numBlocks>{});

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Matrix::N()>{}, [&](const auto i)
        {
            using LSTraits = LinearSolverTraits<i>;
            using DiagBlock = DiagBlockType<Matrix, i>;
            using RHSBlock = VecBlockType<Vector, i>;

            using Smoother = Dune::SeqSSOR<DiagBlock, RHSBlock, RHSBlock>;
            using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
            SmootherArgs args;
            args.iterations = 1;
            args.relaxationFactor = 1;

            Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
            params.setDebugLevel(this->verbosity());
            params.setDefaultValuesIsotropic(LSTraits::GridView::dimension);

            using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<DiagBlock, Dune::Amg::FirstDiagonal>>;
            Criterion criterion(params);

            auto& linearOperator = std::get<i>(linearOperators);
            using SeqLinearOperator = typename LSTraits::template Sequential<DiagBlock, RHSBlock>::LinearOperator;
            linearOperator = std::make_shared<SeqLinearOperator>(m[Dune::index_constant<i>{}][Dune::index_constant<i>{}]);

            // Cast the linear operator from a pointer to the base class
            // to a pointer to the actually employed derived class.
            auto lop = std::dynamic_pointer_cast<SeqLinearOperator>(linearOperator);

            using AMG = Dune::Amg::AMG<SeqLinearOperator, RHSBlock, Smoother>;
            std::get<i>(preconditioners) = std::make_shared<AMG>(*lop, criterion, args);
        });

        TuplePreconditioner<Vector, decltype(preconditioners)> prec(preconditioners);

        Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector> op(m);

        static const int restartGMRes = getParamFromGroup<double>(this->paramGroup(), "LinearSolver.GMResRestart", 10);
        Dune::RestartedGMResSolver<Vector> solver(op, prec, this->residReduction(), restartGMRes,
                                                  this->maxIter(), this->verbosity());

        auto bTmp(b);
        solver.apply(x, bTmp, result_);
    }

    template<template<class X, std::size_t i> class LinearOperator, class Vector, std::size_t... Is>
    auto makeLinearOperators_(std::index_sequence<Is...>)
    {
        return std::make_tuple(std::shared_ptr<LinearOperator<Vector, Is>>()...);
    }

    template<template<class X, std::size_t i> class ScalarProduct, class Vector, std::size_t... Is>
    auto makeScalarProducts_(std::index_sequence<Is...>)
    {
        return std::make_tuple(std::shared_ptr<ScalarProduct<Vector, Is>>()...);
    }

    template<template<std::size_t i> class Comm, std::size_t... Is>
    auto makeComms_(std::index_sequence<Is...>)
    {
        return std::make_tuple(std::shared_ptr<Comm<Is>>()...);
    }

    template<template<class X, std::size_t i> class Preconditioner, class Vector, std::size_t... Is>
    auto makePreconditioners_(std::index_sequence<Is...>)
    {
        return std::make_tuple(std::shared_ptr<Preconditioner<Vector, Is>>()...);
    }

    template<class GridView, class DofMapper, class String, std::size_t... Is>
    BlockDiagAMGGMResSolver(const GridView& gridView,
                               const DofMapper& dofMapper,
                               const String& paramGroup,
                               std::index_sequence<Is...> is)
    : parHelpers_(std::make_tuple(std::make_shared<ParallelHelper<Is>>(std::get<Is>(gridView), std::get<Is>(dofMapper))...))
    , isParallel_(Dune::MPIHelper::getCollectiveCommunication().size() > 1)
    {
        reset();
    }

    ParallelHelperTuple parHelpers_;
    Dune::InverseOperatorResult result_;
    bool isParallel_;
    bool firstCall_;
};

// \}

} // end namespace Dumux

#endif
