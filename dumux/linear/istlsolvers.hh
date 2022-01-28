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
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/linear/algebratraits.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/solvercategory.hh>
#include <dumux/linear/solver.hh>

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

} // end namespace Dumux::Detail


namespace Dumux {

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
    static constexpr bool convertMultiTypeVectorAndMatrix
        = convertMultiTypeLATypes && isMultiTypeBlockVector<XVector>::value;
#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using ScalarProduct = Dune::ScalarProduct<typename InverseOperator::domain_type>;
    using ParallelHelper = ParallelISTLHelper<LinearSolverTraits>;
#endif

public:
    /*!
     * \brief Constructor for sequential solvers
     */
    IstlIterativeLinearSolver(const std::string& paramGroup = "")
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
    template <class GridView, class DofMapper>
    IstlIterativeLinearSolver(const GridView& gridView,
                              const DofMapper& dofMapper,
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
    template <class GridView, class DofMapper>
    IstlIterativeLinearSolver(std::shared_ptr<Comm> communication,
                     std::shared_ptr<ScalarProduct> scalarProduct,
                     const GridView& gridView,
                     const DofMapper& dofMapper,
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
        if constexpr (convertMultiTypeVectorAndMatrix)
        {
            // create the bcrs matrix the IterativeSolver backend can handle
            auto M = MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A);

            // get the new matrix sizes
            const std::size_t numRows = M.N();
            assert(numRows == M.M());

            // create the vector the IterativeSolver backend can handle
            auto bTmp = VectorConverter<BVector>::multiTypeToBlockVector(b);
            assert(bTmp.size() == numRows);

            // create a blockvector to which the linear solver writes the solution
            using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
            using BlockVector = typename Dune::BlockVector<VectorBlock>;
            BlockVector y(numRows);

            auto linearOperator = std::make_shared<Dune::MatrixAdapter<decltype(M), decltype(y), decltype(bTmp)>>(M);
            auto solver = constructPreconditionedSolver_(linearOperator);

            // solve linear system
            solver.apply(y, bTmp, result_);

            // copy back the result y into x
            if(result_.converged)
                VectorConverter<XVector>::retrieveValues(x, y);
        }
        else
        {
            // construct solver from linear operator
            using SequentialTraits = typename LinearSolverTraits::template Sequential<Matrix, XVector>;
            auto linearOperator = std::make_shared<typename SequentialTraits::LinearOperator>(A);
            auto solver = constructPreconditionedSolver_(linearOperator);

            // solve linear system
            solver.apply(x, b, result_);
        }

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
        using TL = Dune::TypeList<typename LinearOperator::matrix_type, typename LinearOperator::domain_type, typename LinearOperator::range_type>;
        std::shared_ptr<Prec> prec = PreconditionerFactory{}(TL{}, op, params);

#if HAVE_MPI && DUNE_VERSION_GT(DUNE_ISTL,2,7)
        if (prec->category() != op->category() && prec->category() == Dune::SolverCategory::sequential)
            prec = Dune::wrapPreconditioner4Parallel(prec, op);
#endif
        return {op, scalarProduct_, prec, params_};
    }

#if HAVE_MPI
    std::unique_ptr<ParallelHelper> parallelHelper_;
    std::shared_ptr<Comm> communication_;
#endif
    Dune::SolverCategory::Category solverCategory_;
    std::shared_ptr<ScalarProduct> scalarProduct_;

    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::string name_;
};


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
    IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlDefaultBlockLevelPreconditionerFactory<Dune::SeqILU>,
        /*accepts multi-type istl types?*/ true
    >;

/*!
 * \ingroup Linear
 * \brief An AMG preconditioned BiCGSTAB solver using dune-istl
 */
template<class LSTraits, class LATraits>
using AMGBiCGSTABIstlSolver =
    IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::BiCGSTABSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlAmgPreconditionerFactory
    >;

/*!
 * \ingroup Linear
 * \brief An AMG preconditioned CG solver using dune-istl
 */
template<class LSTraits, class LATraits>
using AMGCGIstlSolver =
    IstlIterativeLinearSolver<LSTraits, LATraits,
        Dune::CGSolver<typename LATraits::SingleTypeVector>,
        Detail::IstlAmgPreconditionerFactory
    >;


/*!
 * \ingroup Linear
 * \brief Direct dune-istl linear solvers
 */
template<template<class M> class Solver>
class DirectIstlSolver : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        // support dune-istl multi-type block vector/matrix
        if constexpr (isMultiTypeBlockVector<Vector>())
        {
            auto [AA, xx, bb] = convertIstlMultiTypeToBCRSSystem_(A, x, b);
            bool converged = solve_(AA, xx, bb);
            if (converged)
                VectorConverter<Vector>::retrieveValues(x, xx);
            return converged;
        }

        /*
           Copy vectors into standard istl block vector.
           This is necessary for all model _not_ using a FieldVector<Scalar, blockSize>
           Could be avoided for vectors that already have the right type using another if constexpr branch
           but it shouldn't impact performance too much
         */
        else
        {
            constexpr auto blockSize = Detail::blockSize<decltype(b[0])>();
            using Scalar = std::decay_t<decltype(b[0][0])>;
            using BlockType = Dune::FieldVector<Scalar, blockSize>;
            Dune::BlockVector<BlockType> xTmp; xTmp.resize(b.size());
            Dune::BlockVector<BlockType> bTmp(xTmp);

            bTmp = b;
            bool converged = solve_(A, xTmp, bTmp);
            if (converged)
                x  = xTmp;
            return converged;
        }
    }

    std::string name() const
    {
        return "Direct solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    Dune::InverseOperatorResult result_;

    template<class Matrix, class Vector>
    bool solve_(const Matrix& A, Vector& x, const Vector& b)
    {
        static_assert(isBCRSMatrix<Matrix>::value, "Direct solver only works with BCRS matrices!");
        using BlockType = typename Matrix::block_type;
        static_assert(BlockType::rows == BlockType::cols, "Matrix block must be quadratic!");
        constexpr auto blockSize = BlockType::rows;

        Solver<Matrix> solver(A, this->verbosity() > 0);

        Vector bTmp(b);
        solver.apply(x, bTmp, result_);

        int size = x.size();
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                using std::isnan;
                using std::isinf;
                if (isnan(x[i][j]) || isinf(x[i][j]))
                {
                    result_.converged = false;
                    break;
                }
            }
        }

        return result_.converged;
    }

    template<class Matrix, class Vector>
    auto convertIstlMultiTypeToBCRSSystem_(const Matrix& A, Vector& x, const Vector& b)
    {
        const auto AA = MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A);

        // get the new matrix sizes
        const std::size_t numRows = AA.N();
        assert(numRows == AA.M());

        // create the vector the IterativeSolver backend can handle
        const auto bb = VectorConverter<Vector>::multiTypeToBlockVector(b);
        assert(bb.size() == numRows);

        // create a blockvector to which the linear solver writes the solution
        using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
        using BlockVector = typename Dune::BlockVector<VectorBlock>;
        BlockVector xx(numRows);

        return std::make_tuple(std::move(AA), std::move(xx), std::move(bb));
    }
};

} // end namespace Dumux

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
using SuperLUIstlSolver = DirectIstlSolver<Dune::SuperLU>;

} // end namespace Dumux

#endif // HAVE_SUPERLU

#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Direct linear solver using the UMFPack library.
 *
 * See: Davis, Timothy A. (2004). "Algorithm 832". ACM Transactions on
 * Mathematical Software 30 (2): 196–199. doi:10.1145/992200.992206.
 * http://faculty.cse.tamu.edu/davis/suitesparse.html
 */
using UMFPackIstlSolver = DirectIstlSolver<Dune::UMFPack>;

} // end namespace Dumux

#endif // HAVE_UMFPACK

#endif
