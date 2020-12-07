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
 * \brief Linear solvers from dune-istl
 */
#ifndef DUMUX_LINEAR_ISTL_SOLVERS_HH
#define DUMUX_LINEAR_ISTL_SOLVERS_HH

#include <memory>

#include <dune/common/version.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/scalarproducts.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelhelpers.hh>

namespace Dumux::Detail {

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

template<class LinearSolverTraits>
Dune::SolverCategory::Category solverCategory(const typename LinearSolverTraits::GridView& gridView)
{
    if constexpr (LinearSolverTraits::canCommunicate)
    {
        if (gridView.comm().size() <= 1)
            return Dune::SolverCategory::sequential;

        if (LinearSolverTraits::isNonOverlapping(gridView))
            return Dune::SolverCategory::nonoverlapping;
        else
            return Dune::SolverCategory::overlapping;
    }
    else
        return Dune::SolverCategory::sequential;
}

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Base class for linear solvers
 */
template<class LinearSolverTraits, class LinearAlgebraTraits,
         class InverseOperator, class Preconditioner>
class IstlLinearSolver
{
    using Matrix = typename LinearAlgebraTraits::Matrix;
    using XVector = typename InverseOperator::domain_type;
    using BVector = typename InverseOperator::range_type;
    using Scalar = typename InverseOperator::real_type;
#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using ScalarProduct = Dune::ScalarProduct<XVector>;
#endif

public:
    /*!
     * \brief Constructor for sequential solvers
     */
    IstlLinearSolver(const std::string& paramGroup = "")
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
            DUNE_THROW(Dune::InvalidStateException, "Using sequential constructor for parallel run. Use signature with gridView and dofMapper!");

        initializeParameters_(paramGroup);
        solverCategory_ = Dune::SolverCategory::sequential;
        scalarProduct_ = std::make_shared<ScalarProduct>();
    }

    /*!
     * \brief Constructor for parallel and sequential solvers
     */
    IstlLinearSolver(const typename LinearSolverTraits::GridView& gridView,
                     const typename LinearSolverTraits::DofMapper& dofMapper,
                     const std::string& paramGroup = "")
    {
        initializeParameters_(paramGroup);
#if HAVE_MPI
        solverCategory_ = Detail::solverCategory<LinearSolverTraits>(gridView);
        if (solverCategory_ != Dune::SolverCategory::sequential)
        {
            parallelHelper_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
            communication_ = std::make_shared<Comm>(gridView.comm(), solverCategory_);
            scalarProduct_ = Dune::createScalarProduct<XVector>(*communication_, solverCategory_);
            parallelHelper_->createParallelIndexSet(*communication_);
        }
        else
            scalarProduct_ = std::make_shared<ScalarProduct>();
#else
        solverCategory_ = Dune::SolverCategory::sequential;
#endif
    }

#if HAVE_MPI
    /*!
     * \brief Constructor with custom scalar product and communication
     */
    IstlLinearSolver(std::shared_ptr<Comm> communication,
                     std::shared_ptr<ScalarProduct> scalarProduct,
                     const typename LinearSolverTraits::GridView& gridView,
                     const typename LinearSolverTraits::DofMapper& dofMapper,
                     const std::string& paramGroup = "")
    {
        initializeParameters_(paramGroup);
        solverCategory_ = Detail::solverCategory(gridView);
        scalarProduct_ = scalarProduct;
        communication_ = communication;
        if (solverCategory_ != Dune::SolverCategory::sequential)
        {
            parallelHelper_ = std::make_unique<ParallelISTLHelper<LinearSolverTraits>>(gridView, dofMapper);
            parallelHelper_->createParallelIndexSet(communication);
        }
    }
#endif

    bool solve(Matrix& A, XVector& x, BVector& b)
    {
#if HAVE_MPI
        return solveSequentialOrParallel_(A, x, b);
#else
        return solveSequential_(A, x, b);
#endif
    }

    Scalar norm(const XVector& x) const
    {
#if HAVE_MPI
        if (solverCategory_ == Dune::SolverCategory::nonoverlapping)
        {
            auto y(x); // make a copy because the vector needs to be made consistent
            using GV = typename LinearSolverTraits::GridView;
            using DM = typename LinearSolverTraits::DofMapper;
            ParallelVectorHelper<GV, DM, LinearSolverTraits::dofCodim> vectorHelper(parallelHelper_->gridView(), parallelHelper_->dofMapper());
            vectorHelper.makeNonOverlappingConsistent(y);
            return scalarProduct_->norm(y);
        }
        else
#endif
        {
            return scalarProduct_->norm(x);
        }
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

    const std::string& name() const
    {
        return name_;
    }

    void setResidualReduction(double residReduction)
    { params_["reduction"] = std::to_string(residReduction); }

private:
    //! reset the linear solver factory
    void initializeParameters_(const std::string& paramGroup)
    {
        params_ = LinearSolverParameters<LinearSolverTraits>::createParameterTree(paramGroup);
    }

    bool solveSequential_(Matrix& A, XVector& x, BVector& b)
    {
        // construct solver from linear operator
        using SequentialTraits = typename LinearSolverTraits::template Sequential<Matrix, XVector>;
        auto linearOperator = std::make_shared<typename SequentialTraits::LinearOperator>(A);
        auto solver = constructPreconditionedSolver_(linearOperator);

        // solve linear system
        solver.apply(x, b, result_);
        return result_.converged;
    }

#if HAVE_MPI
    bool solveSequentialOrParallel_(Matrix& A, XVector& x, BVector& b)
    {
        // Dune::MultiTypeBlockMatrix does not provide a ColIterator which is needed by Dune::NonoverlappingSchwarzOperator.
        // We therefore can only solve these types of systems sequentially.
        // TODO: This can be adapted once the situation in Dune ISTL changes.
        if constexpr (isMultiTypeBlockMatrix<Matrix>::value)
            return solveSequential_(A, x, b);

        else
        {
            switch (solverCategory_)
            {
                case Dune::SolverCategory::sequential:
                    return solveSequential_(A, x, b);
                case Dune::SolverCategory::nonoverlapping:
                    using NOTraits = typename LinearSolverTraits::template ParallelNonoverlapping<Matrix, XVector>;
                    return solveParallel_<NOTraits>(A, x, b);
                case Dune::SolverCategory::overlapping:
                    using OTraits = typename LinearSolverTraits::template ParallelOverlapping<Matrix, XVector>;
                    return solveParallel_<OTraits>(A, x, b);
                default: DUNE_THROW(Dune::InvalidStateException, "Unknown solver category");
            }
        }
    }

    template<class ParallelTraits>
    bool solveParallel_(Matrix& A, XVector& x, BVector& b)
    {
#if DUNE_VERSION_GT(DUNE_ISTL,2,7)
        // make linear algebra consistent
        prepareLinearAlgebraParallel<LinearSolverTraits, ParallelTraits>(A, b, *parallelHelper_);

        // construct solver from linear operator
        auto linearOperator = std::make_shared<typename ParallelTraits::LinearOperator>(A, *communication_);
        auto solver = constructPreconditionedSolver_(linearOperator);

        // solve linear system
        solver.apply(x, b, result_);
        return result_.converged;
#else
        DUNE_THROW(Dune::NotImplemented, "Parallel solvers only available for dune-istl > 2.7");
#endif
    }
#endif // HAVE_MPI

    template<class LinearOperator>
    InverseOperator constructPreconditionedSolver_(std::shared_ptr<LinearOperator>& op)
    {
        const auto& params = params_.sub("preconditioner");
        using Prec = Dune::Preconditioner<typename LinearOperator::domain_type, typename LinearOperator::range_type>;
        std::shared_ptr<Prec> prec = std::make_shared<Preconditioner>(op, params);

#if HAVE_MPI && DUNE_VERSION_GT(DUNE_ISTL,2,7)
        if (prec->category() != op->category() && prec->category() == Dune::SolverCategory::sequential)
            prec = Dune::wrapPreconditioner4Parallel(prec, op);
#endif
        return {op, scalarProduct_, prec, params_};
    }

#if HAVE_MPI
    std::unique_ptr<ParallelISTLHelper<LinearSolverTraits>> parallelHelper_;
    std::shared_ptr<Comm> communication_;
#endif
    Dune::SolverCategory::Category solverCategory_;
    std::shared_ptr<ScalarProduct> scalarProduct_;

    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::string name_;
};

template<class LSTraits, class LATraits>
using ILUBiCGSTABIstlSolver = IstlLinearSolver<LSTraits, LATraits,
                                               Dune::BiCGSTABSolver<typename LATraits::Vector>,
                                               Dune::SeqILU<typename LATraits::Matrix,
                                                            typename LATraits::Vector,
                                                            typename LATraits::Vector,
                                                            Detail::preconditionerBlockLevel<typename LATraits::Matrix>()>
                                               >;

} // end namespace Dumux

#endif
