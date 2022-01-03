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
#include <dumux/linear/solver.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/operators.hh>
#include <dumux/linear/scalarproducts.hh>
#include <dumux/linear/preconditioners.hh>

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
template<class LinearSolverTraits,
         class InverseOperator, class Preconditioner,
         class Matrix = typename Preconditioner::matrix_type>
class IstlLinearSolver
{
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

template<class LinearSolverTraits, class Matrix, class Vector>
using ILUBiCGSTABIstlSolver = IstlLinearSolver<LinearSolverTraits,
                                               Dune::BiCGSTABSolver<Vector>,
                                               Dune::SeqILU<Matrix, Vector, Vector, Detail::preconditionerBlockLevel<Matrix>()>
                                               >;

/*!
 * \ingroup Linear
 * \brief Linear solvers preconditioned by a block-diagonal AMG
 */
template<class LinearSolverTraitsTuple, class InverseOperator, class Matrix>
class BlockDiagAMGPreconditionedSolver : public LinearSolver
{
    using Vector = typename InverseOperator::domain_type;
    using Scalar = typename InverseOperator::real_type;
    using Category = Dune::SolverCategory::Category;
    static constexpr auto numBlocks = std::tuple_size_v<LinearSolverTraitsTuple>;

    template<std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<Vector>()[Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using LinearSolverTraits = std::tuple_element_t<i, LinearSolverTraitsTuple>;

    template<std::size_t i>
    using ScalarProduct = std::shared_ptr<Dune::ScalarProduct<VecBlockType<i>>>;
    using ScalarProductTuple = typename makeFromIndexedType<std::tuple,
                                                            ScalarProduct,
                                                            std::make_index_sequence<numBlocks>
                                                           >::type;

    template<std::size_t i>
    using Comm = std::shared_ptr<Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>>;
    using CommTuple = typename makeFromIndexedType<std::tuple,
                                                   Comm,
                                                   std::make_index_sequence<numBlocks>
                                                  >::type;

    template<std::size_t i>
    using ParallelHelper = ParallelISTLHelper<LinearSolverTraits<i>>;
    template<std::size_t i>
    using ParallelHelperSP = std::shared_ptr<ParallelHelper<i>>;
    using ParallelHelperTuple = typename makeFromIndexedType<std::tuple,
                                                             ParallelHelperSP,
                                                             std::make_index_sequence<numBlocks>
                                                            >::type;

public:
    /*!
     * \brief Construct the linear solver
     *
     * \param gridViews a tuple of grid views
     * \param dofMappers a tuple of dof mappers
     * \param paramGroup parameter group
     */
    template<class GridViewTuple, class DofMapperTuple>
    BlockDiagAMGPreconditionedSolver(const GridViewTuple& gridViews,
                                     const DofMapperTuple& dofMappers,
                                     const std::string& paramGroup = "")
    : BlockDiagAMGPreconditionedSolver(gridViews, dofMappers, paramGroup, std::make_index_sequence<numBlocks>{})
    {}

    /*!
     * \brief Solve a linear system
     *
     * \param A the matrix
     * \param x the seeked solution vector, containing an initial guess
     * \param b the right hand side vector
     */
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        using Prec = BlockDiagAMGPreconditioner<LinearSolverTraitsTuple, Matrix, Vector>;
        auto prec = std::make_shared<Prec>(A, b, comms_, parHelpers_);

        using LOP = TupleLinearOperator<Vector, Matrix, decltype(prec->linearOperators())>;
        auto op = std::make_shared<LOP>(prec->linearOperators(), A);

        auto rank = Dune::MPIHelper::getCollectiveCommunication().rank();
        if (rank != 0)
            params_["verbose"] = "0";
        InverseOperator solver(op, scalarProduct_, prec, params_);

        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    /*!
     * \brief Result containing the convergence history
     */
    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    /*!
     * \brief Calculate the norm of a right-hand side vector
     *
     * \param x vector for which the norm will be calculated
     *
     * The norm calculation employs the scalar product after
     * possibly making the parts consistent that correspond
     * to a non-overlapping distribution in the parallel regime.
     */
    Scalar norm(const Vector& x) const
    {
        auto y(x); // make a copy because the vector needs to be made consistent
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numBlocks>{}, [&](const auto i)
        {
            if (categories_[i] == Dune::SolverCategory::nonoverlapping)
            {
                using GV = typename LinearSolverTraits<i>::GridView;
                using DM = typename LinearSolverTraits<i>::DofMapper;
                using PVHelper = ParallelVectorHelper<GV, DM, LinearSolverTraits<i>::dofCodim>;

                const auto& parHelper = *std::get<i>(parHelpers_);

                PVHelper vectorHelper(parHelper.gridView(), parHelper.dofMapper());

                vectorHelper.makeNonOverlappingConsistent(y[i]);
            }
        });

        return scalarProduct_->norm(y);
    }

    /*!
     * \brief Name of the solver
     */
    std::string name() const
    { return "block-diagonal AMG preconditioned GMRes solver"; }

private:

    template <class ParallelTraits, class Comm, class SP, class PH>
    void prepareCommAndScalarProduct_(std::shared_ptr<Comm>& comm, SP& scalarProduct, PH& parHelper, Category& category)
    {
        if constexpr (ParallelTraits::isNonOverlapping)
            category = Dune::SolverCategory::nonoverlapping;
        else
            category = Dune::SolverCategory::overlapping;

        comm = std::make_shared<Comm>(parHelper.gridView().comm(), category);
        parHelper.createParallelIndexSet(*comm);
        scalarProduct = std::make_shared<typename ParallelTraits::ScalarProduct>(*comm);
    }

    template<class GridViewTuple, class DofMapperTuple, std::size_t... Is>
    BlockDiagAMGPreconditionedSolver(const GridViewTuple& gridViews,
                                     const DofMapperTuple& dofMappers,
                                     const std::string& paramGroup,
                                     std::index_sequence<Is...> is)
    : parHelpers_(std::make_tuple(std::make_shared<ParallelHelper<Is>>(std::get<Is>(gridViews), std::get<Is>(dofMappers))...))
    {
        params_ = LinearSolverParameters<LinearSolverTraits<0>>::createParameterTree(paramGroup);

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<numBlocks>{}, [&](const auto i)
        {
            using LSTraits = LinearSolverTraits<i>;
            using DiagBlock = std::decay_t<decltype(std::declval<Matrix>()[i][i])>;
            using RHSBlock = VecBlockType<i>;

            auto& comm = std::get<i>(comms_);
            auto& scalarProduct = std::get<i>(scalarProducts_);
            auto& parHelper = *std::get<i>(parHelpers_);

            if (LSTraits::isNonOverlapping(parHelper.gridView()))
            {
                using PTraits = typename LSTraits::template ParallelNonoverlapping<DiagBlock, RHSBlock>;
                prepareCommAndScalarProduct_<PTraits>(comm, scalarProduct, parHelper, categories_[i]);
            }
            else
            {
                using PTraits = typename LSTraits::template ParallelOverlapping<DiagBlock, RHSBlock>;
                prepareCommAndScalarProduct_<PTraits>(comm, scalarProduct, parHelper, categories_[i]);
            }
        });

        scalarProduct_ = std::make_shared<TupleScalarProduct<Vector, ScalarProductTuple>>(scalarProducts_);
    }

    ParallelHelperTuple parHelpers_;
    CommTuple comms_;
    ScalarProductTuple scalarProducts_;
    std::shared_ptr<TupleScalarProduct<Vector, ScalarProductTuple>> scalarProduct_;
    std::array<Category, numBlocks> categories_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
};

template<class LinearSolverTraitsTuple, class Matrix, class Vector>
using BlockDiagAMGBiCGSTABSolver = BlockDiagAMGPreconditionedSolver<LinearSolverTraitsTuple,
                                                                    Dune::BiCGSTABSolver<Vector>,
                                                                    Matrix
                                                                   >;

template<class LinearSolverTraitsTuple, class Matrix, class Vector>
using BlockDiagAMGGMResSolver = BlockDiagAMGPreconditionedSolver<LinearSolverTraitsTuple,
                                                                 Dune::RestartedGMResSolver<Vector>,
                                                                 Matrix
                                                                >;

} // end namespace Dumux

#endif
