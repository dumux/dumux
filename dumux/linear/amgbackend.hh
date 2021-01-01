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
 * \brief Provides a parallel linear solver based on the ISTL AMG preconditioner
 *        and the ISTL BiCGSTAB solver.
 */
#ifndef DUMUX_PARALLEL_AMGBACKEND_HH
#define DUMUX_PARALLEL_AMGBACKEND_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvers.hh>

#include <dumux/linear/solver.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/solvercategory.hh>
#include <dumux/linear/istlsolvers.hh>

namespace Dumux {

namespace Detail {

template <class LinearSolverTraits>
class OldAMGBiCGSTABBackend : public LinearSolver
{
public:
    /*!
     * \brief Construct the backend for the sequential case only
     *
     * \param paramGroup the parameter group for parameter lookup
     */
    [[deprecated("Use new AMGBiCGSTABBackend<LinearSolverTraits, LinearAlgebraTraits> with 2nd template parameter.")]]
    OldAMGBiCGSTABBackend(const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , isParallel_(Dune::MPIHelper::getCommunication().size() > 1)
    {
        if (isParallel_)
            DUNE_THROW(Dune::InvalidStateException, "Using sequential constructor for parallel run. Use signature with gridView and dofMapper!");

        checkAvailabilityOfDirectSolver_();
    }

    /*!
     * \brief Construct the backend for parallel or sequential runs
     *
     * \param gridView the grid view on which we are performing the multi-grid
     * \param dofMapper an index mapper for dof entities
     * \param paramGroup the parameter group for parameter lookup
     */
    [[deprecated("Use new AMGBiCGSTABBackend<LinearSolverTraits, LinearAlgebraTraits> with 2nd template parameter.")]]
    OldAMGBiCGSTABBackend(const typename LinearSolverTraits::GridView& gridView,
                          const typename LinearSolverTraits::DofMapper& dofMapper,
                          const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
#if HAVE_MPI
    , isParallel_(Dune::MPIHelper::getCommunication().size() > 1)
#endif
    {
#if HAVE_MPI
        if constexpr (LinearSolverTraits::canCommunicate)
        {
            if (isParallel_)
                phelper_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
        }
#endif
        checkAvailabilityOfDirectSolver_();
    }

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
            phelper_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
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
        return result_.converged;
    }

    /*!
     * \brief The name of the solver
     */
    std::string name() const
    {
        return "AMG-preconditioned BiCGSTAB solver";
    }

    /*!
     * \brief The result containing the convergence history.
     */
    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

    template<class Vector>
    Scalar norm(const Vector& x) const = delete;

private:
    //! see https://gitlab.dune-project.org/core/dune-istl/-/issues/62
    void checkAvailabilityOfDirectSolver_()
    {
#if !HAVE_SUPERLU && !HAVE_UMFPACK
        std::cout << "\nAMGBiCGSTABBackend: No direct solver backend found. Using iterative solver as coarse grid solver.\n"
                  << "Note that dune-istl currently hard-codes a tolerance of 1e-2 for the iterative coarse grid solver.\n"
                  << "This may result in reduced accuracy or performance depending on your setup.\nConsider installing "
                  << "UMFPack (SuiteSparse) or SuperLU or apply the istl patch, see dumux/patches/README.md." << std::endl;
#endif
    }

#if HAVE_MPI
    template<class Matrix, class Vector>
    void solveSequentialOrParallel_(Matrix& A, Vector& x, Vector& b)
    {
        if constexpr (LinearSolverTraits::canCommunicate)
        {
            if (isParallel_)
            {
                if (LinearSolverTraits::isNonOverlapping(phelper_->gridView()))
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

        std::shared_ptr<Comm> comm;
        std::shared_ptr<LinearOperator> linearOperator;
        std::shared_ptr<ScalarProduct> scalarProduct;
        prepareLinearAlgebraParallel<LinearSolverTraits, ParallelTraits>(A, b, comm, linearOperator, scalarProduct, *phelper_);

        using SeqSmoother = Dune::SeqSSOR<Matrix, Vector, Vector>;
        using Smoother = typename ParallelTraits::template Preconditioner<SeqSmoother>;
        solveWithAmg_<Smoother>(A, x, b, linearOperator, comm, scalarProduct);
    }
#endif // HAVE_MPI

    template<class Matrix, class Vector>
    void solveSequential_(Matrix& A, Vector& x, Vector& b)
    {
        using Comm = Dune::Amg::SequentialInformation;
        using Traits = typename LinearSolverTraits::template Sequential<Matrix, Vector>;
        using LinearOperator = typename Traits::LinearOperator;
        using ScalarProduct = typename Traits::ScalarProduct;

        auto comm = std::make_shared<Comm>();
        auto linearOperator = std::make_shared<LinearOperator>(A);
        auto scalarProduct = std::make_shared<ScalarProduct>();

        using Smoother = Dune::SeqSSOR<Matrix, Vector, Vector>;
        solveWithAmg_<Smoother>(A, x, b, linearOperator, comm, scalarProduct);
    }

    template<class Smoother, class Matrix, class Vector, class LinearOperator, class Comm, class ScalarProduct>
    void solveWithAmg_(Matrix& A, Vector& x, Vector& b,
                       std::shared_ptr<LinearOperator>& linearOperator,
                       std::shared_ptr<Comm>& comm,
                       std::shared_ptr<ScalarProduct>& scalarProduct)
    {
        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<Matrix, Dune::Amg::FirstDiagonal>>;

        //! \todo Check whether the default accumulation mode atOnceAccu is needed.
        //! \todo make parameters changeable at runtime from input file / parameter tree
        Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(LinearSolverTraits::GridView::dimension);
        params.setDebugLevel(this->verbosity());
        Criterion criterion(params);
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        using Amg = Dune::Amg::AMG<LinearOperator, Vector, Smoother, Comm>;
        auto amg = std::make_shared<Amg>(*linearOperator, criterion, smootherArgs, *comm);

        Dune::BiCGSTABSolver<Vector> solver(*linearOperator, *scalarProduct, *amg, this->residReduction(), this->maxIter(),
                                            comm->communicator().rank() == 0 ? this->verbosity() : 0);

        solver.apply(x, b, result_);
    }

#if HAVE_MPI
    std::unique_ptr<ParallelISTLHelper<LinearSolverTraits>> phelper_;
#endif
    Dune::InverseOperatorResult result_;
    bool isParallel_ = false;
};

/*!
 * \ingroup Linear
 * \brief An AMG preconditioned BiCGSTAB solver using dune-istl
 */
template<class LSTraits, class LATraits>
using NewAMGBiCGSTABBackend =
    IstlLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlAmgPreconditionerFactory
    >;

template<int i>
struct AMGImplHelper
{
    template<class LSTraits, class LATraits>
    using type = NewAMGBiCGSTABBackend<LSTraits, LATraits>;
};

template<>
struct AMGImplHelper<1>
{
    template<class T>
    using type = OldAMGBiCGSTABBackend<T>;
};

} // end namespace Detail

/*!
 * \ingroup Linear
 * \brief A linear solver based on the ISTL AMG preconditioner
 *        and the ISTL BiCGSTAB solver.
 */
template<class ...T>
using AMGBiCGSTABBackend = typename Detail::AMGImplHelper<sizeof...(T)>::template type<T...>;

} // end namespace Dumux

#endif
