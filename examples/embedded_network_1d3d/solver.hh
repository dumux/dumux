// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Simple ilu0 block diagonal preconditioner that can store decomposition
 */
#ifndef DUMUX_EXAMPLE_1D3D_MULTIDOMAIN_LINEAR_SOLVER_HH
#define DUMUX_EXAMPLE_1D3D_MULTIDOMAIN_LINEAR_SOLVER_HH

#include <type_traits>
#include <tuple>
#include <utility>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelmatrixadapter.hh>

namespace Dumux::Example {

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
template<class MatrixType, class VectorType>
class BlockDiagILU0BiCGSTABSolver : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<class Matrix, class Vector>
    bool solve(const Matrix& m, Vector& x, const Vector& b)
    {
        setup(m);
        return solve(x, b);
    }

    void setup(const MatrixType& m)
    {
        preconditioner_ = std::make_unique<BlockDiagILU0Preconditioner<MatrixType, VectorType, VectorType>>(m);
        linearOperator_ = std::make_shared<Dumux::ParallelMultiTypeMatrixAdapter<MatrixType, VectorType, VectorType>>(m);
        solver_ = std::make_unique<Dune::BiCGSTABSolver<VectorType>>(*linearOperator_, *preconditioner_, this->residReduction(), this->maxIter(), this->verbosity());

        isSetup_ = true;
    }

    bool solve(VectorType& x, const VectorType& b)
    {
        assert(isSetup_);

        auto bTmp(b);
        solver_->apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    std::string name() const
    { return "block-diagonal ILU0-preconditioned BiCGSTAB solver"; }

private:
    bool isSetup_ = false;

    std::shared_ptr<Dune::MatrixAdapter<MatrixType, VectorType, VectorType>> linearOperator_;
    std::unique_ptr<Dune::Preconditioner<VectorType, VectorType>> preconditioner_;
    std::unique_ptr<Dune::InverseOperator<VectorType, VectorType>> solver_;

    Dune::InverseOperatorResult result_;
};

} // end namespace Dumux::Example

#endif
