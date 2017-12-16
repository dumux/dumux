// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Schur complement solver for the embedded problem, the idea is to solve
 *        only linear systems of the size of the much smaller embedded domain by
 *        forming a Schur complement preconditioner. If the embedded system is much
 *        smaller it might be feasible to just use a direct solver.
 */
#ifndef DUMUX_EMBEDDED_SCHUR_COMPLEMENT_SOLVER_BACKEND_HH
#define DUMUX_EMBEDDED_SCHUR_COMPLEMENT_SOLVER_BACKEND_HH

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/basicproperties.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dumux/linear/amgproperties.hh>

namespace Dumux {

/*!
  \brief Schur complement linear operator
  this applies the Schur complement to a vector and can be used in an inverse operator
*/
template<class AType, class BType, class CType, class DType, class InvOpDType, class X, class Y>
class SchurComplement : public Dune::LinearOperator<X,Y>
{
public:
    // export types
    // using matrix_type = DType;
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;

    //! constructor: just store a reference to a matrix
    explicit SchurComplement (const AType& A, const BType& B,
                              const CType& C, const DType& D,
                              const std::shared_ptr<InvOpDType>& invD)
    : A_(A), B_(B), C_(C), D_(D), invD_(invD) {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
        // apply C, note x and aTmp1 have different size (Cx)
        X aTmp1(D_.N());
        aTmp1 = 0.0;
        C_.mv(x, aTmp1);

        // apply D^-1 (D^-1Cx)
        auto aTmp2 = aTmp1;
        Dune::InverseOperatorResult result;
        invD_->apply(aTmp2, aTmp1, result);

        // apply B (BD^-1Cx)
        B_.mv(aTmp2, y);

        y *= -1.0;

        // add Ax (Ax - BD^-1Cx) -> the full Schur complement of S_D := D/M
        A_.umv(x, y);
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
        auto tmp = y;
        this->apply(x, tmp);
        y.axpy(alpha, tmp);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
      return Dune::SolverCategory::sequential;
    }

  private:
    const AType& A_;
    const BType& B_;
    const CType& C_;
    const DType& D_;
    std::shared_ptr<InvOpDType> invD_;
};

/*! \brief The schur complement preconditioner

   \tparam M The matrix type to operate on
   \tparam X Type of the update
   \tparam Y Type of the defect
   \tparam l The block level to invert. Default is 1
*/
template<class TypeTag, class AType, class BType, class CType, class DType, class InvOpDType, class X, class Y>
class SchurComplementPreconditioner : public Dune::Preconditioner<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    // using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    using BulkVector = typename GET_PROP(TypeTag, SolutionVector)::SolutionVectorBulk;
    using EmbeddedVector = typename GET_PROP(TypeTag, SolutionVector)::SolutionVectorLowDim;

    /*! \brief Approximates the Schur complement
     */
    SchurComplementPreconditioner (const AType& A, const BType& B,
                                   const CType& C, const DType& D,
                                   const std::shared_ptr<InvOpDType>& invD)
    : A_(A), B_(B), C_(C), D_(D), invD_(invD) {}


    virtual void pre (X& x, Y& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
       \brief Apply the preconditioner.

       \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply (X& x, const Y& b)
    {
        using namespace Dune::Indices;
        const BulkVector& f = b[_0];
        const EmbeddedVector& g = b[_1];

        static const int verbosity = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, PreconditionerVerbosity);
        static const int schurMaxIter = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, SchurIterations);
        static const double schurReduction = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, SchurResidualReduction);

        Dune::Richardson<BulkVector, BulkVector> identity(1.0);
        SchurComplement<AType, BType, CType, DType, InvOpDType, BulkVector, BulkVector> schurOp(A_, B_, C_, D_, invD_);
        Dune::BiCGSTABSolver<BulkVector> invOpS(schurOp, identity, schurReduction, schurMaxIter, verbosity);

        // apply S^-1 to solve the bulk block
        BulkVector fTmp(f);
        Dune::InverseOperatorResult res;
        invOpS.apply(x[_0], fTmp, res);

        // prepare rhs for application of D^-1 (g - Ca)
        EmbeddedVector gTmp(g);
        C_.mmv(x[_0], gTmp);

        // then solve for embedded block b using the action of D^-1
        // use direct solver
        invD_->apply(x[_1], gTmp, res);
    }

    virtual void post (X& x)
    { DUNE_UNUSED_PARAMETER(x); }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    { return Dune::SolverCategory::sequential; }

private:
    const AType& A_;
    const BType& B_;
    const CType& C_;
    const DType& D_;
    std::shared_ptr<InvOpDType> invD_;
};

template <class TypeTag>
class SchurComplementSolver
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

public:
    SchurComplementSolver() = default;
    SchurComplementSolver(const Problem& problem) {}

    // expects a system as a multi-type matrix
    // | A  B |
    // | C  D |

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<class Matrix, class Vector>
    bool solve(const Matrix& M, Vector& x, const Vector& b)
    {
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        static const double reduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        static const double maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        static const int restartGMRes = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, GMResRestart);

        using AType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockBulk;
        using BType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockBulkCoupling;
        using CType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockLowDimCoupling;
        using DType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockLowDim;

        using namespace Dune::Indices;
        const AType& A = M[_0][_0];
        const BType& B = M[_0][_1];
        const CType& C = M[_1][_0];
        const DType& D = M[_1][_1];

        // Solver to compute applications of D^-1 used in the preconditioner
        // We construct it here so e.g. LU decomposition is only done once
        using InnerSolver = Dune::UMFPack<DType>;
        auto invD = std::make_shared<InnerSolver>(D, false);

        SchurComplementPreconditioner<TypeTag, AType, BType, CType, DType, InnerSolver, Vector, Vector> preconditioner(A, B, C, D, invD);

        Dune::MatrixAdapter<Matrix, Vector, Vector> op(M);
        Dune::RestartedfGMResSolver<Vector> solver(op, preconditioner, reduction, restartGMRes, maxIter, verbosity);
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return true;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};


template<class FirstBlock, class FirstX, class FirstY,
         class SecondBlock, class SecondX, class SecondY,
         class X, class Y>
class BlockILU0Preconditioner : public Dune::Preconditioner<X, Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    // using matrix_type = typename std::remove_const<M>::type;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param A The matrix to operate on.
       \param w The relaxation factor.
     */
    BlockILU0Preconditioner (const FirstBlock& A, const SecondBlock& B)
      : firstILU_(A, 1.0), secondILU_(B, 1.0) // does ILU decomposition
    {}

    /*!
       \brief Prepare the preconditioner.

       \copydoc Preconditioner::pre(X&,Y&)
     */
    virtual void pre (X& x, Y& b)
    {
      DUNE_UNUSED_PARAMETER(x);
      DUNE_UNUSED_PARAMETER(b);
    }

    /*!
       \brief Apply the preconditoner.

       \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply (X& v, const Y& d)
    {
        using namespace Dune::Indices;
        firstILU_.apply(v[_0], d[_0]);
        secondILU_.apply(v[_1], d[_1]);
    }

    /*!
       \brief Clean up.

       \copydoc Preconditioner::post(X&)
     */
    virtual void post (X& x)
    {
      DUNE_UNUSED_PARAMETER(x);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
      return Dune::SolverCategory::sequential;
    }

  private:
    Dune::SeqILU0<FirstBlock, FirstX, FirstY> firstILU_;
    Dune::SeqILU0<SecondBlock, SecondX, SecondY> secondILU_;
};

template <class TypeTag>
class BlockILU0BiCGSTABSolver
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

public:
    BlockILU0BiCGSTABSolver() = default;
    BlockILU0BiCGSTABSolver(const Problem& problem) {}

    // expects a system as a multi-type matrix
    // | A  B |
    // | C  D |

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<class Matrix, class Vector>
    bool solve(const Matrix& M, Vector& x, const Vector& b)
    {
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        static const double reduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        static const double maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);

        using AType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockBulk;
        using DType = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockLowDim;

        using BulkVector = typename GET_PROP(TypeTag, SolutionVector)::SolutionVectorBulk;
        using EmbeddedVector = typename GET_PROP(TypeTag, SolutionVector)::SolutionVectorLowDim;

        using namespace Dune::Indices;
        const AType& A = M[_0][_0];
        const DType& D = M[_1][_1];

        // Solver to compute applications of D^-1 used in the preconditioner
        // We construct it here so e.g. LU decomposition is only done once
        using InnerSolver = Dune::UMFPack<DType>;
        auto invD = std::make_shared<InnerSolver>(D, false);

        BlockILU0Preconditioner<AType, BulkVector, BulkVector, DType, EmbeddedVector, EmbeddedVector, Vector, Vector> preconditioner(A, D);
        Dune::MatrixAdapter<Matrix, Vector, Vector> op(M);
        Dune::BiCGSTABSolver<Vector> solver(op, preconditioner, reduction, maxIter, verbosity);
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};

} // end namespace Dumux

#endif
