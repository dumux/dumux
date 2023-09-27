// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Linear solvers from dune-istl
 */
#ifndef DUMUX_LINEAR_ISTL_SOLVERS_HH
#define DUMUX_LINEAR_ISTL_SOLVERS_HH

#include <memory>
#include <variant>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/solverfactory.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/solvercategory.hh>
#include <dumux/linear/solver.hh>

#include <dune/istl/foreach.hh>

namespace Dumux::Detail::IstlSolvers {

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

template<template<class,class,class,int> class Preconditioner, int blockLevel = 1>
class IstlDefaultBlockLevelPreconditionerFactory
{
public:
    template<class TL, class M>
    auto operator() (TL typeList, const M& matrix, const Dune::ParameterTree& config)
    {
        using Matrix = typename Dune::TypeListElement<0, decltype(typeList)>::type;
        using Domain = typename Dune::TypeListElement<1, decltype(typeList)>::type;
        using Range = typename Dune::TypeListElement<2, decltype(typeList)>::type;
        std::shared_ptr<Dune::Preconditioner<Domain, Range>> preconditioner
            = std::make_shared<Preconditioner<Matrix, Domain, Range, blockLevel>>(matrix, config);
        return preconditioner;
    }
};

template<template<class,class,class> class Preconditioner>
class IstlDefaultPreconditionerFactory
{
    template<class TL, class M>
    auto operator() (TL typeList, const M& matrix, const Dune::ParameterTree& config)
    {
        using Matrix = typename Dune::TypeListElement<0, decltype(typeList)>::type;
        using Domain = typename Dune::TypeListElement<1, decltype(typeList)>::type;
        using Range = typename Dune::TypeListElement<2, decltype(typeList)>::type;
        std::shared_ptr<Dune::Preconditioner<Domain, Range>> preconditioner
            = std::make_shared<Preconditioner<Matrix, Domain, Range>>(matrix, config);
        return preconditioner;
    }
};

using IstlAmgPreconditionerFactory = Dune::AMGCreator;

template<class M, bool convert = false>
struct MatrixForSolver { using type = M; };

template<class M>
struct MatrixForSolver<M, true>
{ using type = std::decay_t<decltype(MatrixConverter<M>::multiTypeToBCRSMatrix(std::declval<M>()))>; };

template<class V, bool convert = false>
struct VectorForSolver { using type = V; };

template<class V>
struct VectorForSolver<V, true>
{ using type = std::decay_t<decltype(VectorConverter<V>::multiTypeToBlockVector(std::declval<V>()))>; };

template<class LSTraits, class LATraits, bool convert, bool parallel = LSTraits::canCommunicate>
struct MatrixOperator;

template<class LSTraits, class LATraits, bool convert>
struct MatrixOperator<LSTraits, LATraits, convert, true>
{
    using M = typename MatrixForSolver<typename LATraits::Matrix, convert>::type;
    using V = typename VectorForSolver<typename LATraits::Vector, convert>::type;
#if HAVE_MPI
    using type = std::variant<
        std::shared_ptr<typename LSTraits::template Sequential<M, V>::LinearOperator>,
        std::shared_ptr<typename LSTraits::template ParallelOverlapping<M, V>::LinearOperator>,
        std::shared_ptr<typename LSTraits::template ParallelNonoverlapping<M, V>::LinearOperator>
    >;
#else
    using type = std::variant<
        std::shared_ptr<typename LSTraits::template Sequential<M, V>::LinearOperator>
    >;
#endif
};

template<class LSTraits, class LATraits, bool convert>
struct MatrixOperator<LSTraits, LATraits, convert, false>
{
    using M = typename MatrixForSolver<typename LATraits::Matrix, convert>::type;
    using V = typename VectorForSolver<typename LATraits::Vector, convert>::type;
    using type = std::variant<
        std::shared_ptr<typename LSTraits::template Sequential<M, V>::LinearOperator>
    >;
};

} // end namespace Dumux::Detail::IstlSolvers

namespace Dumux::Detail {

struct IstlSolverResult : public Dune::InverseOperatorResult
{
    IstlSolverResult() = default;
    IstlSolverResult(const IstlSolverResult&) = default;
    IstlSolverResult(IstlSolverResult&&) = default;

    IstlSolverResult(const Dune::InverseOperatorResult& o) : InverseOperatorResult(o) {}
    IstlSolverResult(Dune::InverseOperatorResult&& o) : InverseOperatorResult(std::move(o)) {}

    operator bool() const { return this->converged; }
};

/*!
 * \ingroup Linear
 * \brief Standard dune-istl iterative linear solvers
 */
template<class LinearSolverTraits, class LinearAlgebraTraits,
         class InverseOperator, class PreconditionerFactory,
         bool convertMultiTypeLATypes = false>
class IstlIterativeLinearSolver
{
    using Matrix = typename LinearAlgebraTraits::Matrix;
    using XVector = typename LinearAlgebraTraits::Vector;
    using BVector = typename LinearAlgebraTraits::Vector;
    using Scalar = typename InverseOperator::real_type;

    using ScalarProduct = Dune::ScalarProduct<typename InverseOperator::domain_type>;

    static constexpr bool convertMultiTypeVectorAndMatrix
        = convertMultiTypeLATypes && isMultiTypeBlockVector<XVector>::value;
    using MatrixForSolver = typename Detail::IstlSolvers::MatrixForSolver<Matrix, convertMultiTypeVectorAndMatrix>::type;
    using BVectorForSolver = typename Detail::IstlSolvers::VectorForSolver<BVector, convertMultiTypeVectorAndMatrix>::type;
    using XVectorForSolver = typename Detail::IstlSolvers::VectorForSolver<XVector, convertMultiTypeVectorAndMatrix>::type;
    // a variant type that can hold sequential, overlapping, and non-overlapping operators
    using MatrixOperatorHolder = typename Detail::IstlSolvers::MatrixOperator<
        LinearSolverTraits, LinearAlgebraTraits, convertMultiTypeVectorAndMatrix
    >::type;

#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using ParallelHelper = ParallelISTLHelper<LinearSolverTraits>;
#endif

    using ParameterInitializer = std::variant<std::string, Dune::ParameterTree>;
public:

    /*!
     * \brief Constructor for sequential solvers
     */
    IstlIterativeLinearSolver(const ParameterInitializer& params = "")
    {
        if (Dune::MPIHelper::getCommunication().size() > 1)
            DUNE_THROW(Dune::InvalidStateException, "Using sequential constructor for parallel run. Use signature with gridView and dofMapper!");

        initializeParameters_(params);
        solverCategory_ = Dune::SolverCategory::sequential;
        scalarProduct_ = std::make_shared<ScalarProduct>();
    }

    /*!
     * \brief Constructor for parallel and sequential solvers
     */
    template <class GridView, class DofMapper>
    IstlIterativeLinearSolver(const GridView& gridView,
                              const DofMapper& dofMapper,
                              const ParameterInitializer& params = "")
    {
        initializeParameters_(params);
#if HAVE_MPI
        solverCategory_ = Detail::solverCategory<LinearSolverTraits>(gridView);
        if constexpr (LinearSolverTraits::canCommunicate)
        {

            if (solverCategory_ != Dune::SolverCategory::sequential)
            {
                parallelHelper_ = std::make_shared<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
                communication_ = std::make_shared<Comm>(gridView.comm(), solverCategory_);
                scalarProduct_ = Dune::createScalarProduct<XVector>(*communication_, solverCategory_);
                parallelHelper_->createParallelIndexSet(*communication_);
            }
            else
                scalarProduct_ = std::make_shared<ScalarProduct>();
        }
        else
            scalarProduct_ = std::make_shared<ScalarProduct>();
#else
        solverCategory_ = Dune::SolverCategory::sequential;
        scalarProduct_ = std::make_shared<ScalarProduct>();
#endif
    }

#if HAVE_MPI
    /*!
     * \brief Constructor with custom scalar product and communication
     */
    template <class GridView, class DofMapper>
    IstlIterativeLinearSolver(std::shared_ptr<Comm> communication,
                              std::shared_ptr<ScalarProduct> scalarProduct,
                              const GridView& gridView,
                              const DofMapper& dofMapper,
                              const ParameterInitializer& params = "")
    {
        initializeParameters_(params);
        solverCategory_ = Detail::solverCategory(gridView);
        scalarProduct_ = scalarProduct;
        communication_ = communication;
        if constexpr (LinearSolverTraits::canCommunicate)
        {
            if (solverCategory_ != Dune::SolverCategory::sequential)
            {
                parallelHelper_ = std::make_shared<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
                parallelHelper_->createParallelIndexSet(communication);
            }
        }
    }
#endif

    /*!
     * \brief Solve the linear system Ax = b
     */
    IstlSolverResult solve(Matrix& A, XVector& x, BVector& b)
    { return solveSequentialOrParallel_(A, x, b); }

    /*!
     * \brief Set the matrix A of the linear system Ax = b for reuse
     */
    void setMatrix(std::shared_ptr<Matrix> A)
    {
        linearOperator_ = makeParallelOrSequentialLinearOperator_(std::move(A));
        solver_ = constructPreconditionedSolver_(linearOperator_);
    }

    /*!
     * \brief Set the matrix A of the linear system Ax = b for reuse
     * \note The client has to take care of the lifetime management of A
     */
    void setMatrix(Matrix& A)
    { setMatrix(Dune::stackobject_to_shared_ptr(A)); }

    /*!
     * \brief Solve the linear system Ax = b where A has been set with \ref setMatrix
     */
    IstlSolverResult solve(XVector& x, BVector& b) const
    {
        if (!solver_)
            DUNE_THROW(Dune::InvalidStateException, "Called solve(x, b) but no linear operator has been set");

        return solveSequentialOrParallel_(x, b, *solver_);
    }

    /*!
     * \brief Compute the 2-norm of vector x
     */
    Scalar norm(const XVector& x) const
    {
#if HAVE_MPI
        if constexpr (LinearSolverTraits::canCommunicate)
        {
            if (solverCategory_ == Dune::SolverCategory::nonoverlapping)
            {
                auto y(x); // make a copy because the vector needs to be made consistent
                using GV = typename LinearSolverTraits::GridView;
                using DM = typename LinearSolverTraits::DofMapper;
                ParallelVectorHelper<GV, DM, LinearSolverTraits::dofCodim> vectorHelper(parallelHelper_->gridView(), parallelHelper_->dofMapper());
                vectorHelper.makeNonOverlappingConsistent(y);
                return scalarProduct_->norm(y);
            }
        }
#endif
        if constexpr (convertMultiTypeVectorAndMatrix)
        {
            auto y = VectorConverter<XVector>::multiTypeToBlockVector(x);
            return scalarProduct_->norm(y);
        }
        else
            return scalarProduct_->norm(x);
    }

    /*!
     * \brief The name of the linear solver
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Set the residual reduction tolerance
     */
    void setResidualReduction(double residReduction)
    {
        params_["reduction"] = std::to_string(residReduction);

        // reconstruct the solver with new parameters
        if (solver_)
            solver_ = constructPreconditionedSolver_(linearOperator_);
    }

    /*!
     * \brief Set the maximum number of linear solver iterations
     */
    void setMaxIter(std::size_t maxIter)
    {
        params_["maxit"] = std::to_string(maxIter);

        // reconstruct the solver with new parameters
        if (solver_)
            solver_ = constructPreconditionedSolver_(linearOperator_);
    }

    /*!
     * \brief Set the linear solver parameters
     * \param params Either a std::string giving a parameter group (parameters are read from input file) or a Dune::ParameterTree
     * \note In case of a Dune::ParameterTree, the parameters are passed trait to the linear solver and preconditioner
     */
    void setParams(const ParameterInitializer& params)
    {
        initializeParameters_(params);

        // reconstruct the solver with new parameters
        if (solver_)
            solver_ = constructPreconditionedSolver_(linearOperator_);
    }

private:

    void initializeParameters_(const ParameterInitializer& params)
    {
        if (std::holds_alternative<std::string>(params))
            params_ = Dumux::LinearSolverParameters<LinearSolverTraits>::createParameterTree(std::get<std::string>(params));
        else
            params_ = std::get<Dune::ParameterTree>(params);
    }

    MatrixOperatorHolder makeSequentialLinearOperator_(std::shared_ptr<Matrix> A)
    {
        using SequentialTraits = typename LinearSolverTraits::template Sequential<MatrixForSolver, XVectorForSolver>;
        if constexpr (convertMultiTypeVectorAndMatrix)
        {
            // create the BCRS matrix the IterativeSolver backend can handle
            auto M = std::make_shared<MatrixForSolver>(MatrixConverter<Matrix>::multiTypeToBCRSMatrix(*A));
            return std::make_shared<typename SequentialTraits::LinearOperator>(M);
        }
        else
        {
            return std::make_shared<typename SequentialTraits::LinearOperator>(A);
        }
    }

    template<class ParallelTraits>
    MatrixOperatorHolder makeParallelLinearOperator_(std::shared_ptr<Matrix> A, ParallelTraits = {})
    {
#if HAVE_MPI
        // make matrix consistent
        prepareMatrixParallel<LinearSolverTraits, ParallelTraits>(*A, *parallelHelper_);
        return std::make_shared<typename ParallelTraits::LinearOperator>(std::move(A), *communication_);
#else
        DUNE_THROW(Dune::InvalidStateException, "Calling makeParallelLinearOperator for sequential run");
#endif
    }

    MatrixOperatorHolder makeParallelOrSequentialLinearOperator_(std::shared_ptr<Matrix> A)
    {
        return executeSequentialOrParallel_(
            [&]{ return makeSequentialLinearOperator_(std::move(A)); },
            [&](auto traits){ return makeParallelLinearOperator_(std::move(A), traits); }
        );
    }

    MatrixOperatorHolder makeSequentialLinearOperator_(Matrix& A)
    { return makeSequentialLinearOperator_(Dune::stackobject_to_shared_ptr<Matrix>(A)); }

    MatrixOperatorHolder makeParallelOrSequentialLinearOperator_(Matrix& A)
    { return makeParallelOrSequentialLinearOperator_(Dune::stackobject_to_shared_ptr<Matrix>(A)); }

    template<class ParallelTraits>
    MatrixOperatorHolder makeParallelLinearOperator_(Matrix& A, ParallelTraits = {})
    { return makeParallelLinearOperator_<ParallelTraits>(Dune::stackobject_to_shared_ptr<Matrix>(A)); }

    IstlSolverResult solveSequential_(Matrix& A, XVector& x, BVector& b)
    {
        // construct solver from linear operator
        auto linearOperatorHolder = makeSequentialLinearOperator_(A);
        auto solver = constructPreconditionedSolver_(linearOperatorHolder);

        return solveSequential_(x, b, *solver);
    }

    IstlSolverResult solveSequential_(XVector& x, BVector& b, InverseOperator& solver) const
    {
        Dune::InverseOperatorResult result;
        if constexpr (convertMultiTypeVectorAndMatrix)
        {
            // create the vector the IterativeSolver backend can handle
            BVectorForSolver bTmp = VectorConverter<BVector>::multiTypeToBlockVector(b);

            // create a block vector to which the linear solver writes the solution
            XVectorForSolver y(bTmp.size());

            // solve linear system
            solver.apply(y, bTmp, result);

            // copy back the result y into x
            if (result.converged)
                VectorConverter<XVector>::retrieveValues(x, y);
        }
        else
        {
            // solve linear system
            solver.apply(x, b, result);
        }

        return result;
    }

    IstlSolverResult solveSequentialOrParallel_(Matrix& A, XVector& x, BVector& b)
    {
        return executeSequentialOrParallel_(
            [&]{ return solveSequential_(A, x, b); },
            [&](auto traits){ return solveParallel_(A, x, b, traits); }
        );
    }

    IstlSolverResult solveSequentialOrParallel_(XVector& x, BVector& b, InverseOperator& solver) const
    {
        return executeSequentialOrParallel_(
            [&]{ return solveSequential_(x, b, solver); },
            [&](auto traits){ return solveParallel_(x, b, solver, traits); }
        );
    }

    template<class ParallelTraits>
    IstlSolverResult solveParallel_(Matrix& A, XVector& x, BVector& b, ParallelTraits = {})
    {
        // construct solver from linear operator
        auto linearOperatorHolder = makeParallelLinearOperator_<ParallelTraits>(A);
        auto solver = constructPreconditionedSolver_(linearOperatorHolder);
        return solveParallel_<ParallelTraits>(x, b, *solver);
    }

    template<class ParallelTraits>
    IstlSolverResult solveParallel_(XVector& x, BVector& b, InverseOperator& solver, ParallelTraits = {}) const
    {
#if HAVE_MPI
        // make right hand side consistent
        prepareVectorParallel<LinearSolverTraits, ParallelTraits>(b, *parallelHelper_);

        // solve linear system
        Dune::InverseOperatorResult result;
        solver.apply(x, b, result);
        return result;
#else
        DUNE_THROW(Dune::InvalidStateException, "Calling makeParallelLinearOperator for sequential run");
#endif
    }


    std::shared_ptr<InverseOperator> constructPreconditionedSolver_(MatrixOperatorHolder& ops)
    {
        return std::visit([&](auto&& op)
        {
            using LinearOperator = typename std::decay_t<decltype(op)>::element_type;
            const auto& params = params_.sub("preconditioner");
            using Prec = Dune::Preconditioner<typename LinearOperator::domain_type, typename LinearOperator::range_type>;
            using TL = Dune::TypeList<typename LinearOperator::matrix_type, typename LinearOperator::domain_type, typename LinearOperator::range_type>;
            std::shared_ptr<Prec> prec = PreconditionerFactory{}(TL{}, op, params);

#if HAVE_MPI
            if (prec->category() != op->category() && prec->category() == Dune::SolverCategory::sequential)
                prec = Dune::wrapPreconditioner4Parallel(prec, op);
#endif
            return std::make_shared<InverseOperator>(op, scalarProduct_, prec, params_);
        }, ops);
    }

    template<class Seq, class Par>
    decltype(auto) executeSequentialOrParallel_(Seq&& sequentialAction, Par&& parallelAction) const
    {
#if HAVE_MPI
        // For Dune::MultiTypeBlockMatrix there is currently no generic way
        // of handling parallelism, we therefore can only solve these types of systems sequentially
        if constexpr (isMultiTypeBlockMatrix<Matrix>::value || !LinearSolverTraits::canCommunicate)
            return sequentialAction();
        else
        {
            switch (solverCategory_)
            {
                case Dune::SolverCategory::sequential:
                    return sequentialAction();
                case Dune::SolverCategory::nonoverlapping:
                    using NOTraits = typename LinearSolverTraits::template ParallelNonoverlapping<Matrix, XVector>;
                    return parallelAction(NOTraits{});
                case Dune::SolverCategory::overlapping:
                    using OTraits = typename LinearSolverTraits::template ParallelOverlapping<Matrix, XVector>;
                    return parallelAction(OTraits{});
                default: DUNE_THROW(Dune::InvalidStateException, "Unknown solver category");
            }
        }
#else
        return sequentialAction();
#endif
    }

#if HAVE_MPI
    std::shared_ptr<const ParallelHelper> parallelHelper_;
    std::shared_ptr<Comm> communication_;
#endif

    Dune::SolverCategory::Category solverCategory_;
    std::shared_ptr<ScalarProduct> scalarProduct_;

    // for stored solvers (reuse matrix)
    MatrixOperatorHolder linearOperator_;
    // for stored solvers (reuse matrix)
    std::shared_ptr<InverseOperator> solver_;

    Dune::ParameterTree params_;
    std::string name_;
};

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief An ILU preconditioned BiCGSTAB solver using dune-istl
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. It can be applied to
 * nonsymmetric matrices.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n indicates
 * fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template<class LSTraits, class LATraits>
using ILUBiCGSTABIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlSolvers::IstlDefaultBlockLevelPreconditionerFactory<Dune::SeqILU>,
        // the Dune::ILU preconditioners don't accept multi-type matrices
        /*convert multi-type istl types?*/ true
    >;

/*!
 * \ingroup Linear
 * \brief An ILU preconditioned GMres solver using dune-istl
 *
 * Solver: The GMRes (generalized minimal residual) method is an iterative
 * method for the numerical solution of a nonsymmetric system of linear
 * equations.\n
 * See: Saad, Y., Schultz, M. H. (1986). "GMRES: A generalized minimal residual
 * algorithm for solving nonsymmetric linear systems." SIAM J. Sci. and Stat.
 * Comput. 7: 856–869.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n indicates
 * fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template<class LSTraits, class LATraits>
using ILURestartedGMResIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::RestartedGMResSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlSolvers::IstlDefaultBlockLevelPreconditionerFactory<Dune::SeqILU>,
        // the Dune::ILU preconditioners don't accept multi-type matrices
        /*convert multi-type istl types?*/ true
    >;

/*!
 * \ingroup Linear
 * \brief An SSOR-preconditioned BiCGSTAB solver using dune-istl
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template<class LSTraits, class LATraits>
using SSORBiCGSTABIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::Vector>,
        Detail::IstlSolvers::IstlDefaultBlockLevelPreconditionerFactory<Dune::SeqSSOR>
    >;

/*!
 * \ingroup Linear
 * \brief An SSOR-preconditioned CG solver using dune-istl
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
template<class LSTraits, class LATraits>
using SSORCGIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::CGSolver<typename LATraits::Vector>,
        Detail::IstlSolvers::IstlDefaultBlockLevelPreconditionerFactory<Dune::SeqSSOR>
    >;

/*!
 * \ingroup Linear
 * \brief An AMG preconditioned BiCGSTAB solver using dune-istl
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: AMG (algebraic multigrid)
 */
template<class LSTraits, class LATraits>
using AMGBiCGSTABIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlSolvers::IstlAmgPreconditionerFactory,
        // the AMG preconditioner doesn't accept multi-type matrices
        /*convert multi-type istl types?*/ true
    >;

/*!
 * \ingroup Linear
 * \brief An AMG preconditioned CG solver using dune-istl
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: AMG (algebraic multigrid)
 */
template<class LSTraits, class LATraits>
using AMGCGIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::CGSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlSolvers::IstlAmgPreconditionerFactory,
        // the AMG preconditioner doesn't accept multi-type matrices
        /*convert multi-type istl types?*/ true
    >;

/*!
 * \ingroup Linear
 * \brief An Uzawa preconditioned BiCGSTAB solver using dune-istl
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: Uzawa method for saddle point problems
 * \note Expects a 2x2 MultiTypeBlockMatrix
 */
template<class LSTraits, class LATraits>
using UzawaBiCGSTABIstlSolver =
    Detail::IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::Vector>,
        Detail::IstlSolvers::IstlDefaultBlockLevelPreconditionerFactory<Dumux::SeqUzawa>
    >;

} // end namespace Dumux

namespace Dumux::Detail {

/*!
 * \ingroup Linear
 * \brief Direct dune-istl linear solvers
 */
template<class LSTraits, class LATraits, template<class M> class Solver,
         bool convertMultiTypeVectorAndMatrix = isMultiTypeBlockVector<typename LATraits::Vector>::value>
class DirectIstlSolver : public LinearSolver
{
    using Matrix = typename LATraits::Matrix;
    using XVector = typename LATraits::Vector;
    using BVector = typename LATraits::Vector;

    using MatrixForSolver = typename Detail::IstlSolvers::MatrixForSolver<Matrix, convertMultiTypeVectorAndMatrix>::type;
    using BVectorForSolver = typename Detail::IstlSolvers::VectorForSolver<BVector, convertMultiTypeVectorAndMatrix>::type;
    using XVectorForSolver = typename Detail::IstlSolvers::VectorForSolver<XVector, convertMultiTypeVectorAndMatrix>::type;
    using InverseOperator = Dune::InverseOperator<XVectorForSolver, BVectorForSolver>;
public:
    using LinearSolver::LinearSolver;

    /*!
     * \brief Solve the linear system Ax = b
     */
    IstlSolverResult solve(const Matrix& A, XVector& x, const BVector& b)
    {
        return solve_(A, x, b);
    }

    /*!
     * \brief Solve the linear system Ax = b using the matrix set with \ref setMatrix
     */
    IstlSolverResult solve(XVector& x, const BVector& b)
    {
        if (!solver_)
            DUNE_THROW(Dune::InvalidStateException, "Called solve(x, b) but no linear operator has been set");

        return solve_(x, b, *solver_);
    }

    /*!
     * \brief Set the matrix A of the linear system Ax = b for reuse
     */
    void setMatrix(std::shared_ptr<Matrix> A)
    {
        if constexpr (convertMultiTypeVectorAndMatrix)
            matrix_ = std::make_shared<MatrixForSolver>(MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A));
        else
            matrix_ = A;

        solver_ = std::make_shared<Solver<MatrixForSolver>>(*matrix_);
    }

    /*!
     * \brief Set the matrix A of the linear system Ax = b for reuse
     * \note The client has to take care of the lifetime management of A
     */
    void setMatrix(Matrix& A)
    { setMatrix(Dune::stackobject_to_shared_ptr(A)); }

    /*!
     * \brief name of the linear solver
     */
    std::string name() const
    {
        return "Direct solver";
    }

private:
    IstlSolverResult solve_(const Matrix& A, XVector& x, const BVector& b)
    {
        // support dune-istl multi-type block vector/matrix by copying
        if constexpr (convertMultiTypeVectorAndMatrix)
        {
            const auto AA = MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A);
            Solver<MatrixForSolver> solver(AA, this->verbosity() > 0);
            return solve_(x, b, solver);
        }
        else
        {
            Solver<MatrixForSolver> solver(A, this->verbosity() > 0);
            return solve_(x, b, solver);
        }
    }

    IstlSolverResult solve_(XVector& x, const BVector& b, InverseOperator& solver) const
    {
        Dune::InverseOperatorResult result;

        if constexpr (convertMultiTypeVectorAndMatrix)
        {
            auto bb = VectorConverter<BVector>::multiTypeToBlockVector(b);
            XVectorForSolver xx(bb.size());
            solver.apply(xx, bb, result);
            checkResult_(xx, result);
            if (result.converged)
                VectorConverter<XVector>::retrieveValues(x, xx);
            return result;
        }
        else
        {
            BVectorForSolver bTmp(b);
            solver.apply(x, bTmp, result);
            checkResult_(x, result);
            return result;
        }
    }


    void checkResult_(XVectorForSolver& x, Dune::InverseOperatorResult& result) const
    {
        flatVectorForEach(x, [&](auto&& entry, std::size_t){
            using std::isnan, std::isinf;
            if (isnan(entry) || isinf(entry))
                result.converged = false;
        });
    }

    //! matrix when using the setMatrix interface for matrix reuse
    std::shared_ptr<MatrixForSolver> matrix_;
    //! solver when using the setMatrix interface for matrix reuse
    std::shared_ptr<InverseOperator> solver_;
};

} // end namespace Dumux::Detail

#if HAVE_SUPERLU
#include <dune/istl/superlu.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Direct linear solver using the SuperLU library.
 *
 * See: Li, X. S. (2005). "An overview of SuperLU: Algorithms, implementation,
 * and user interface." ACM Transactions on Mathematical Software (TOMS) 31(3): 302-325.
 * http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
 */
template<class LSTraits, class LATraits>
using SuperLUIstlSolver = Detail::DirectIstlSolver<LSTraits, LATraits, Dune::SuperLU>;

} // end namespace Dumux

#endif // HAVE_SUPERLU

#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#include <dune/common/version.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Direct linear solver using the UMFPack library.
 *
 * See: Davis, Timothy A. (2004). "Algorithm 832". ACM Transactions on
 * Mathematical Software 30 (2): 196–199. doi:10.1145/992200.992206.
 * http://faculty.cse.tamu.edu/davis/suitesparse.html
 */
template<class LSTraits, class LATraits>
using UMFPackIstlSolver = Detail::DirectIstlSolver<
    LSTraits, LATraits, Dune::UMFPack
#if DUNE_VERSION_GTE(DUNE_ISTL,2,10)
    , false // no need to convert multi-type matrix anymore
#endif
>;

} // end namespace Dumux

#endif // HAVE_UMFPACK

#endif
