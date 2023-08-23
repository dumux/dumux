// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Provides a generic linear solver based on the ISTL that chooses the
 *        solver and preconditioner at runtime
 */

#ifndef DUMUX_LINEAR_ISTL_SOLVERFACTORYBACKEND_HH
#define DUMUX_LINEAR_ISTL_SOLVERFACTORYBACKEND_HH

#include <memory>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>

#include "linearsolverparameters.hh"

// preconditioners
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include "preconditioners.hh"

// solvers
#include <dune/istl/solvers.hh>
#include <dune/istl/solverfactory.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/istlsolverregistry.hh>
#include <dumux/linear/solvercategory.hh>

namespace Dumux {

// initSolverFactoriesForMultiTypeBlockMatrix differs in different compilation units,
// so we have it in an anonymous namespace
namespace {

/*!
 * \brief Initializes the direct solvers, preconditioners and iterative solvers in
 *        the factories with the corresponding Matrix and Vector types.
 * \note  We currently consider only direct solvers and preconditioners provided by
 *        Dumux which, unlike the ones implemented in Dune, also support MultiTypeBlockMatrices.
 * \note This function could be removed once Dune::initSolverFactories supports MultiTypeBlockMatrices.
 * \tparam LinearOperator the linear operator type
 */
template<class LinearOperator>
int initSolverFactoriesForMultiTypeBlockMatrix()
{
    using M  = typename LinearOperator::matrix_type;
    using X  = typename LinearOperator::range_type;
    using Y  = typename LinearOperator::domain_type;
    using TL = Dune::TypeList<M,X,Y>;
    auto& dsfac = Dune::DirectSolverFactory<M,X,Y>::instance();
    Dune::addRegistryToFactory<TL>(dsfac, Dumux::MultiTypeBlockMatrixDirectSolverTag{});
    auto& pfac = Dune::PreconditionerFactory<LinearOperator,X,Y>::instance();
    Dune::addRegistryToFactory<TL>(pfac, Dumux::MultiTypeBlockMatrixPreconditionerTag{});
    using TLS = Dune::TypeList<X,Y>;
    auto& isfac = Dune::IterativeSolverFactory<X,Y>::instance();
    return Dune::addRegistryToFactory<TLS>(isfac, Dune::IterativeSolverTag{});
}
} // end namespace

/*!
 * \brief Initialize the solver factories for regular matrices or MultiTypeBlockMatrices
 * \tparam Matrix the matrix
 * \tparam LinearOperator the linear operator
 *
 * \note This function could be removed once Dune::initSolverFactories supports MultiTypeBlockMatrices.
 */
template<class Matrix, class LinearOperator>
void initSolverFactories()
{
    if constexpr (isMultiTypeBlockMatrix<Matrix>::value)
        initSolverFactoriesForMultiTypeBlockMatrix<LinearOperator>();
    else
        Dune::initSolverFactories<LinearOperator>();
}

/*!
 * \ingroup Linear
 * \brief A linear solver using the dune-istl solver factory
 *        to choose the solver and preconditioner at runtime.
 * \note the solvers are configured via the input file
 */
namespace Detail {

template <class LinearSolverTraits>
class OldIstlSolverFactoryBackend : public LinearSolver
{
public:

    /*!
     * \brief Update the solver after grid adaption
     *
     * \param gridView the grid view on which we are performing the multi-grid
     * \param dofMapper an index mapper for dof entities
     */
    void updateAfterGridAdaption(const typename LinearSolverTraits::GridView& gridView,
                                 const typename LinearSolverTraits::DofMapper& dofMapper)
    {
#if HAVE_MPI
        if (isParallel_)
            parallelHelper_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
#endif
    }

    /*!
     * \brief Solve a linear system.
     *
     * \param A the matrix
     * \param x the seeked solution vector, containing the initial solution upon entry
     * \param b the right hand side vector
     */
    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
#if HAVE_MPI
        solveSequentialOrParallel_(A, x, b);
#else
        solveSequential_(A, x, b);
#endif
        firstCall_ = false;
        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

    const std::string& name() const
    {
        return name_;
    }

private:

    void initializeParameters_()
    {
        params_ = Dumux::LinearSolverParameters<LinearSolverTraits>::createParameterTree(paramGroup_);
        checkMandatoryParameters_();
        name_ = params_.get<std::string>("preconditioner.type") + "-preconditioned " + params_.get<std::string>("type");
        if (params_.get<int>("verbose", 0) > 0)
            std::cout << "Initialized linear solver of type: " << name_ << std::endl;
    }

    void checkMandatoryParameters_()
    {
        if (!params_.hasKey("type"))
            DUNE_THROW(Dune::InvalidStateException, "Solver factory needs \"LinearSolver.Type\" parameter to select the solver");

        if (!params_.hasKey("preconditioner.type"))
            DUNE_THROW(Dune::InvalidStateException, "Solver factory needs \"LinearSolver.Preconditioner.Type\" parameter to select the preconditioner");
    }

#if HAVE_MPI
    template<class Matrix, class Vector>
    void solveSequentialOrParallel_(Matrix& A, Vector& x, Vector& b)
    {
        // Dune::MultiTypeBlockMatrix does not provide a ColIterator which is needed by Dune::NonoverlappingSchwarzOperator.
        // We therefore can only solve these types of systems sequentially.
        // TODO: This can be adapted once the situation in Dune ISTL changes.
        if constexpr (LinearSolverTraits::canCommunicate && !isMultiTypeBlockMatrix<Matrix>::value)
        {
            if (isParallel_)
            {
                if (LinearSolverTraits::isNonOverlapping(parallelHelper_->gridView()))
                {
                    using PTraits = typename LinearSolverTraits::template ParallelNonoverlapping<Matrix, Vector>;
                    solveParallel_<PTraits>(A, x, b);
                }
                else
                {
                    using PTraits = typename LinearSolverTraits::template ParallelOverlapping<Matrix, Vector>;
                    solveParallel_<PTraits>(A, x, b);
                }
            }
            else
                solveSequential_(A, x, b);
        }
        else
        {
            solveSequential_(A, x, b);
        }
    }

    template<class ParallelTraits, class Matrix, class Vector>
    void solveParallel_(Matrix& A, Vector& x, Vector& b)
    {
        using Comm = typename ParallelTraits::Comm;
        using LinearOperator = typename ParallelTraits::LinearOperator;
        using ScalarProduct = typename ParallelTraits::ScalarProduct;

        if (firstCall_)
            initSolverFactories<Matrix, LinearOperator>();

        std::shared_ptr<Comm> comm;
        std::shared_ptr<LinearOperator> linearOperator;
        std::shared_ptr<ScalarProduct> scalarProduct;
        prepareLinearAlgebraParallel<LinearSolverTraits, ParallelTraits>(A, b, comm, linearOperator, scalarProduct, *parallelHelper_);

        // construct solver
        auto solver = getSolverFromFactory_(linearOperator);

        // solve linear system
        solver->apply(x, b, result_);
    }
#endif // HAVE_MPI

    template<class Matrix, class Vector>
    void solveSequential_(Matrix& A, Vector& x, Vector& b)
    {
        // construct linear operator
        using Traits = typename LinearSolverTraits::template Sequential<Matrix, Vector>;
        using LinearOperator = typename Traits::LinearOperator;
        auto linearOperator = std::make_shared<LinearOperator>(A);

        if (firstCall_)
            initSolverFactories<Matrix, LinearOperator>();

        // construct solver
        auto solver = getSolverFromFactory_(linearOperator);

        // solve linear system
        solver->apply(x, b, result_);
    }

    template<class LinearOperator>
    auto getSolverFromFactory_(std::shared_ptr<LinearOperator>& fop)
    {
        try { return Dune::getSolverFromFactory(fop, params_); }
        catch(Dune::Exception& e)
        {
            std::cerr << "Could not create solver with factory" << std::endl;
            std::cerr << e.what() << std::endl;
            throw e;
        }
    }

    const std::string paramGroup_;
#if HAVE_MPI
    std::unique_ptr<ParallelISTLHelper<LinearSolverTraits>> parallelHelper_;
#endif
    bool isParallel_ = false;
    bool firstCall_;

    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::string name_;
};

template <class LinearSolverTraits, class LinearAlgebraTraits>
class IstlSolverFactoryBackend : public LinearSolver
{
    using Matrix = typename LinearAlgebraTraits::Matrix;
    using Vector = typename LinearAlgebraTraits::Vector;
    using ScalarProduct = Dune::ScalarProduct<Vector>;
#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
#endif
public:

    /*!
     * \brief Construct the backend for the sequential case only
     *
     * \param paramGroup the parameter group for parameter lookup
     */
    IstlSolverFactoryBackend(const std::string& paramGroup = "")
    : paramGroup_(paramGroup)
    , isParallel_(Dune::MPIHelper::getCommunication().size() > 1)
    {
        if (isParallel_)
            DUNE_THROW(Dune::InvalidStateException, "Using sequential constructor for parallel run. Use signature with gridView and dofMapper!");

        firstCall_ = true;
        initializeParameters_();
        scalarProduct_ = std::make_shared<ScalarProduct>();
        solverCategory_ = Dune::SolverCategory::sequential;
    }

    /*!
     * \brief Construct the backend for parallel or sequential runs
     *
     * \param gridView the grid view for parallel communication via the grid
     * \param dofMapper an index mapper for dof entities
     * \param paramGroup the parameter group for parameter lookup
     */
    IstlSolverFactoryBackend(const typename LinearSolverTraits::GridView& gridView,
                             const typename LinearSolverTraits::DofMapper& dofMapper,
                             const std::string& paramGroup = "")
    : paramGroup_(paramGroup)
#if HAVE_MPI
    , isParallel_(gridView.comm().size() > 1)
#endif
    {
        firstCall_ = true;
        initializeParameters_();
#if HAVE_MPI
        solverCategory_ = Detail::solverCategory<LinearSolverTraits>(gridView);
        if (solverCategory_ != Dune::SolverCategory::sequential)
        {
            parallelHelper_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
            comm_ = std::make_shared<Comm>(gridView.comm(), solverCategory_);
            scalarProduct_ = Dune::createScalarProduct<Vector>(*comm_, solverCategory_);
            parallelHelper_->createParallelIndexSet(*comm_);
        }
        else
            scalarProduct_ = std::make_shared<ScalarProduct>();
#else
        solverCategory_ = Dune::SolverCategory::sequential;
        scalarProduct_ = std::make_shared<ScalarProduct>();
#endif
    }

    /*!
     * \brief Solve a linear system.
     *
     * \param A the matrix
     * \param x the seeked solution vector, containing the initial solution upon entry
     * \param b the right hand side vector
     */
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
#if HAVE_MPI
        solveSequentialOrParallel_(A, x, b);
#else
        solveSequential_(A, x, b);
#endif
        firstCall_ = false;
        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

    const std::string& name() const
    {
        return name_;
    }

    Scalar norm(const Vector& x) const
    {
#if HAVE_MPI
        if constexpr (LinearSolverTraits::canCommunicate && !isMultiTypeBlockVector<Vector>::value)
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
        return scalarProduct_->norm(x);
    }

private:

    void initializeParameters_()
    {
        params_ = Dumux::LinearSolverParameters<LinearSolverTraits>::createParameterTree(paramGroup_);
        checkMandatoryParameters_();
        name_ = params_.get<std::string>("preconditioner.type") + "-preconditioned " + params_.get<std::string>("type");
        if (params_.get<int>("verbose", 0) > 0)
            std::cout << "Initialized linear solver of type: " << name_ << std::endl;
    }

    void checkMandatoryParameters_()
    {
        if (!params_.hasKey("type"))
            DUNE_THROW(Dune::InvalidStateException, "Solver factory needs \"LinearSolver.Type\" parameter to select the solver");

        if (!params_.hasKey("preconditioner.type"))
            DUNE_THROW(Dune::InvalidStateException, "Solver factory needs \"LinearSolver.Preconditioner.Type\" parameter to select the preconditioner");
    }

#if HAVE_MPI
    void solveSequentialOrParallel_(Matrix& A, Vector& x, Vector& b)
    {
        // Dune::MultiTypeBlockMatrix does not provide a ColIterator which is needed by Dune::NonoverlappingSchwarzOperator.
        // We therefore can only solve these types of systems sequentially.
        // TODO: This can be adapted once the situation in Dune ISTL changes.
        if constexpr (LinearSolverTraits::canCommunicate && !isMultiTypeBlockMatrix<Matrix>::value)
        {
            if (isParallel_)
            {
                if (LinearSolverTraits::isNonOverlapping(parallelHelper_->gridView()))
                {
                    using PTraits = typename LinearSolverTraits::template ParallelNonoverlapping<Matrix, Vector>;
                    solveParallel_<PTraits>(A, x, b);
                }
                else
                {
                    using PTraits = typename LinearSolverTraits::template ParallelOverlapping<Matrix, Vector>;
                    solveParallel_<PTraits>(A, x, b);
                }
            }
            else
                solveSequential_(A, x, b);
        }
        else
        {
            solveSequential_(A, x, b);
        }
    }

    template<class ParallelTraits>
    void solveParallel_(Matrix& A, Vector& x, Vector& b)
    {
        using LinearOperator = typename ParallelTraits::LinearOperator;

        if (firstCall_)
            initSolverFactories<Matrix, LinearOperator>();

        prepareLinearAlgebraParallel<LinearSolverTraits, ParallelTraits>(A, b, *parallelHelper_);
        auto linearOperator = std::make_shared<LinearOperator>(A, *comm_);

        // construct solver
        auto solver = getSolverFromFactory_(linearOperator);

        // solve linear system
        solver->apply(x, b, result_);
    }
#endif // HAVE_MPI

    void solveSequential_(Matrix& A, Vector& x, Vector& b)
    {
        // construct linear operator
        using Traits = typename LinearSolverTraits::template Sequential<Matrix, Vector>;
        using LinearOperator = typename Traits::LinearOperator;
        auto linearOperator = std::make_shared<LinearOperator>(A);

        if (firstCall_)
            initSolverFactories<Matrix, LinearOperator>();

        // construct solver
        auto solver = getSolverFromFactory_(linearOperator);

        // solve linear system
        solver->apply(x, b, result_);
    }

    template<class LinearOperator>
    auto getSolverFromFactory_(std::shared_ptr<LinearOperator>& fop)
    {
        try { return Dune::getSolverFromFactory(fop, params_); }
        catch(Dune::Exception& e)
        {
            std::cerr << "Could not create solver with factory" << std::endl;
            std::cerr << e.what() << std::endl;
            throw e;
        }
    }

    const std::string paramGroup_;
#if HAVE_MPI
    std::unique_ptr<ParallelISTLHelper<LinearSolverTraits>> parallelHelper_;
    std::shared_ptr<Comm> comm_;
#endif
    bool isParallel_ = false;
    bool firstCall_;

    std::shared_ptr<ScalarProduct> scalarProduct_;
    Dune::SolverCategory::Category solverCategory_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::string name_;
};

template<int i>
struct IstlSolverFactoryBackendImplHelper
{
    template<class ...T>
    using type = IstlSolverFactoryBackend<T...>;
};

template<>
struct IstlSolverFactoryBackendImplHelper<1>
{
    template<class T>
    using type = OldIstlSolverFactoryBackend<T>;
};

} // end namespace Detail

/*!
 * \ingroup Linear
 * \brief A linear solver using the dune-istl solver factory
 *        to choose the solver and preconditioner at runtime.
 * \note the solvers are configured via the input file
 * \note The template arguments are LinearSolverTraits, LinearAlgebraTraits
 */
template<class ...T>
using IstlSolverFactoryBackend = typename Detail::IstlSolverFactoryBackendImplHelper<sizeof...(T)>::template type<T...>;

} // end namespace Dumux

#endif
