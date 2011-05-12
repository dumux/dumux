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

namespace Dumux {

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
    template <class SeqPreconditioner, class Solver>
    bool solve(const Matrix &M, 
            Vector &x,
            const Vector &b,
            double residReduction,
            int verbosityLevel = 0)
    {
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
        //SeqPreconditioner seqPreCond(*overlapMatrix_, 1, 1.0);
        SeqPreconditioner seqPreCond(*overlapMatrix_, 1.0);
        typedef Dumux::OverlappingPreconditioner<SeqPreconditioner, Overlap> OverlappingPreconditioner;
        OverlappingPreconditioner preCond(seqPreCond, overlapMatrix_->overlap());

        // create the scalar products and linear operators for ISTL
        OverlappingScalarProduct scalarProd(overlapMatrix_->overlap());
        OverlappingOperator opA(*overlapMatrix_);

        // create the actual solver
        Solver solver(opA, 
                scalarProd,
                preCond,
                residReduction,
                /*maxIterations=*/250,
                verbosityLevel);

        // run the solver
        Dune::InverseOperatorResult result;
        solver.apply(*overlapx_, *overlapb_, result);

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
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;

public:
    template <class Problem>
    BoxBiCGStabILU0Solver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b,
            double residReduction,
            int verbosityLevel = 0)
    {
        return ParentType::template solve<SeqPreconditioner, Solver>(M, x, b, residReduction, verbosityLevel);
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
    typedef Dune::CGSolver<OverlappingVector> Solver;

public:
    template <class Problem>
    BoxCGILU0Solver(const Problem &problem, int overlapSize = 3)
    : ParentType(problem, overlapSize)
    {}

    bool solve(const Matrix &M,
            Vector &x,
            const Vector &b,
            double residReduction,
            int verbosityLevel = 0)
    {
        return ParentType::template solve<SeqPreconditioner, Solver>(M, x, b, residReduction, verbosityLevel);
    }
};

} // namespace Dumux

#endif
