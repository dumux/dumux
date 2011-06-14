/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Provides a linear solver for the stabilized BiCG method with
 *        an ILU-0 preconditioner.
 */
#ifndef DUMUX_BOXLINEARSOLVER_HH
#define DUMUX_BOXLINEARSOLVER_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/linear/vertexborderlistfromgrid.hh>
#include <dumux/linear/overlappingbcrsmatrix.hh>
#include <dumux/linear/overlappingblockvector.hh>
#include <dumux/linear/overlappingpreconditioner.hh>
#include <dumux/linear/overlappingscalarproduct.hh>
#include <dumux/linear/overlappingoperator.hh>

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/superlu.hh>

namespace Dumux {

namespace Properties
{
/*!
 * \brief Specifies the verbosity of the linear solver
 *
 * By default it is 0, i.e. it doesn't print anything. Setting this
 * property to 1 prints aggregated convergence rates, 2 prints the
 * convergence rate of every iteration of the scheme.
 */
NEW_PROP_TAG(LSVerbosity);
// the outdated name
NEW_PROP_TAG(NewtonLinearSolverVerbosity);

//! target reduction of the initial residual
// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
NEW_PROP_TAG(LSResidualReduction);

//! maximum number of iterations of solver
NEW_PROP_TAG(LSMaxIterations);

//! relaxation parameter for the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

//! number of preconditioner iterations per solver iteration
NEW_PROP_TAG(PreconditionerIterations);

//! restart parameter for GMRes
NEW_PROP_TAG(GMResRestart);

SET_PROP_DEFAULT(NewtonLinearSolverVerbosity)
{public:
    static constexpr int value = 0;
};

SET_PROP_DEFAULT(LSVerbosity)
{public:
    static constexpr int value = GET_PROP_VALUE(TypeTag, PTAG(NewtonLinearSolverVerbosity));
};

SET_PROP_DEFAULT(LSResidualReduction)
{public:
    static constexpr double value = 1e-6;
};

SET_PROP_DEFAULT(LSMaxIterations)
{public:
    static constexpr int value = 250;
};

SET_PROP_DEFAULT(PreconditionerRelaxation)
{public:
    static constexpr double value = 1.0;
};

SET_PROP_DEFAULT(PreconditionerIterations)
{public:
    static constexpr int value = 1;
};

SET_PROP_DEFAULT(GMResRestart)
{public:
    static constexpr int value = 10;
};
}

/*!
 * \brief Provides a linear solver for the stabilized BiCG method with
 *        an ILU-0 preconditioner.
 *
 * This solver's intention is to be used in conjunction with the box
 * method, so it assumes that the vertices are the only DOFs.
 */
template <class TypeTag>
class BoxLinearSolver
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dumux::OverlappingScalarProduct<OverlappingVector, Overlap> OverlappingScalarProduct;
    typedef Dumux::OverlappingOperator<OverlappingMatrix, OverlappingVector, OverlappingVector> OverlappingOperator;

public:
    BoxLinearSolver(const Problem &problem, int overlapSize)
    : problem_(problem)
    , overlapSize_(overlapSize)
    {
        overlapMatrix_ = 0;
        overlapb_ = 0;
        overlapx_ = 0;
    };

    ~BoxLinearSolver()
    { cleanup_(); }

    /*!
     * \brief Set the structure of the linear system of equations to be solved.
     * 
     * This method allocates space an does the necessary
     * communication before actually calling the solve() method.  As
     * long as the structure of the linear system does not change, the
     * solve method can be called arbitrarily often.
     */
    void setStructureMatrix(const Matrix &M)
    {
        cleanup_();
        prepare_();
    };

    /*!
     * \brief Actually solve the linear system of equations. 
     *
     * \return true if the residual reduction could be achieved, else false.
     */
    template <class PrecBackend, class SolverBackend>
    bool solve(const Matrix &M, 
            Vector &x,
            const Vector &b)
    {
        int verbosity = 0;
        if (problem_.gridView().comm().rank() == 0)
            verbosity = GET_PROP_VALUE(TypeTag, PTAG(LSVerbosity));
        static const int maxIter = GET_PROP_VALUE(TypeTag, PTAG(LSMaxIterations));
        static const double residReduction = GET_PROP_VALUE(TypeTag, PTAG(LSResidualReduction));

        if (!overlapMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        };

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlapMatrix_->assignAdd(M);
        overlapb_->assignAdd(b);
        (*overlapx_) = 0.0;

        // create sequential and overlapping preconditioners
        PrecBackend seqPreCond(*overlapMatrix_);
        typedef typename PrecBackend::Implementation SeqPreconditioner;
        typedef Dumux::OverlappingPreconditioner<SeqPreconditioner, Overlap> OverlappingPreconditioner;
        OverlappingPreconditioner preCond(seqPreCond.imp(), overlapMatrix_->overlap());

        // create the scalar products and linear operators for ISTL
        OverlappingScalarProduct scalarProd(overlapMatrix_->overlap());
        OverlappingOperator opA(*overlapMatrix_);

        // create the actual solver
        SolverBackend solver(opA,
                scalarProd,
                preCond,
                residReduction,
                maxIter,
                verbosity);

        // run the solver
        Dune::InverseOperatorResult result;
        solver.imp().apply(*overlapx_, *overlapb_, result);

        // copy the result back to the non-overlapping vector
        overlapx_->assignTo(x);

        // return the result of the solver
        return result.converged;
    };

private:
    void prepare_(const Matrix &M)
    {
        VertexBorderListFromGrid<GridView, VertexMapper>
        borderListCreator(problem_.gridView(), problem_.vertexMapper());

        // create the overlapping Jacobian matrix
        overlapMatrix_ = new OverlappingMatrix (M,
                borderListCreator.foreignBorderList(),
                borderListCreator.domesticBorderList(),
                overlapSize_);

        // create the overlapping vectors for the residual and the
        // solution
        overlapb_ = new OverlappingVector(overlapMatrix_->overlap());
        overlapx_ = new OverlappingVector(*overlapb_);
    };

    void cleanup_()
    {
        // create the overlapping Jacobian matrix and vectors
        delete overlapMatrix_;
        delete overlapb_;
        delete overlapx_;

        overlapMatrix_ = 0;
        overlapb_ = 0;
        overlapx_ = 0;
    };

    const Problem &problem_;

    int overlapSize_;
    OverlappingMatrix *overlapMatrix_;
    OverlappingVector *overlapb_;
    OverlappingVector *overlapx_;
};

template <class TypeTag, class Imp>
class PrecNoIterBackend
{
public:
    typedef Imp Implementation;

    template <class Matrix>
    PrecNoIterBackend(Matrix& A)
    : imp_(A, GET_PROP_VALUE(TypeTag, PTAG(PreconditionerRelaxation)))
    {}

    Imp& imp()
    {
        return imp_;
    }

private:
    Imp imp_;
};

template <class TypeTag, class Imp>
class PrecIterBackend
{
public:
    typedef Imp Implementation;

    template <class Matrix>
    PrecIterBackend(Matrix& A)
    : imp_(A,
           GET_PROP_VALUE(TypeTag, PTAG(PreconditionerIterations)),
           GET_PROP_VALUE(TypeTag, PTAG(PreconditionerRelaxation)))
    {}

    Imp& imp()
    {
        return imp_;
    }

private:
    Imp imp_;
};

template <class TypeTag, class Imp>
class StandardSolverBackend
{
public:
    typedef Imp Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    template <class Operator, class ScalarProduct, class Prec>
    StandardSolverBackend(Operator& A, ScalarProduct& sp, Prec& prec,
                          Scalar residReduction, int maxIter, int verbosity)
    : imp_(A, sp, prec, residReduction, maxIter, verbosity)
    {}

    Imp& imp()
    {
        return imp_;
    }

private:
    Imp imp_;
};

template <class TypeTag>
class BoxBiCGStabILU0Solver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqILU0<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecNoIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxBiCGStabILU0Solver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxBiCGStabSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxBiCGStabSORSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxBiCGStabSSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxBiCGStabSSORSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxBiCGStabJacSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqJac<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxBiCGStabJacSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxBiCGStabGSSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqGS<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxBiCGStabGSSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxCGILU0Solver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqILU0<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecNoIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxCGILU0Solver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxCGSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxCGSORSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxCGSSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxCGSSORSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxCGJacSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqJac<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxCGJacSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

template <class TypeTag>
class BoxCGGSSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(JacobianMatrix)) Matrix;
    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) Vector;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqGS<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    template <class Problem>
    BoxCGGSSolver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

#if HAVE_SUPERLU
template <class TypeTag>
class SuperLUBackend
{
  typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

public:

  SuperLUBackend(const Problem& problem, int overlapSize = 0)
  : problem_(problem)
  {}

  template<class Matrix, class Vector>
  int solve(const Matrix& A, Vector& x, const Vector& b)
  {
    int verbosity = GET_PROP_VALUE(TypeTag, PTAG(LSVerbosity));
    static const int maxIter = GET_PROP_VALUE(TypeTag, PTAG(LSMaxIterations));
    static const double residReduction = GET_PROP_VALUE(TypeTag, PTAG(LSResidualReduction));

    Vector bTmp(b);

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    enum {numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef typename Dune::BCRSMatrix<MatrixBlock> ISTLMatrix;

    Dune::SuperLU<ISTLMatrix> solver(A);

    solver.apply(x, bTmp, result_);

    return result_.converged;
  }

  const Dune::InverseOperatorResult& result() const
  {
    return result_;
  }

private:
  Dune::InverseOperatorResult result_;
  const Problem& problem_;
};
#endif

} // namespace Dumux

#endif
