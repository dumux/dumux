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
 *
 * \brief Provides a linear solver for the stabilized BiCG method with
 *        an ILU-0 preconditioner.
 */
#ifndef DUMUX_BICGSTAB_ILU0_SOLVER_HH
#define DUMUX_BICGSTAB_ILU0_SOLVER_HH

#include <dumux/common/propertysystem.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/elementborderlistfromgrid.hh>
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/linear/overlappingbcrsmatrix.hh>
#include <dumux/linear/overlappingblockvector.hh>
#include <dumux/linear/overlappingpreconditioner.hh>
#include <dumux/linear/overlappingscalarproduct.hh>
#include <dumux/linear/overlappingoperator.hh>

namespace Dumux {
namespace Properties {
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(PressureCoefficientMatrix);
NEW_PROP_TAG(PressureRHSVector);
NEW_PROP_TAG(SolutionTypes);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(LinearSolverVerbosity);
NEW_PROP_TAG(LinearSolverMaxIterations);
NEW_PROP_TAG(LinearSolverResidualReduction);
NEW_PROP_TAG(PreconditionerRelaxation);
}

/*!
 * \ingroup Linear
 * \brief Provides a linear solver for the stabilized BiCG method with
 *        an ILU-0 preconditioner.
 *
 * This solver's intention is to be used in conjunction with cell centered finite volume (ccvf)
 * methods, so it assumes that the cells are the only DOFs.
 */
template <class TypeTag>
class IMPETBiCGStabILU0Solver
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PressureRHSVector) Vector;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ElementMapper ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Dumux::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef Dumux::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqILU0<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    //typedef Dune::SeqJac<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef Dumux::OverlappingPreconditioner<SeqPreconditioner, Overlap> OverlappingPreconditioner;
    typedef Dumux::OverlappingScalarProduct<OverlappingVector, Overlap> OverlappingScalarProduct;
    typedef Dumux::OverlappingOperator<OverlappingMatrix, OverlappingVector, OverlappingVector> OverlappingOperator;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;

public:
    DUNE_DEPRECATED_MSG("Use ILU0BiCGSTABBackend/AMGBackend from seqsolverbackend.hh/amgbackend.hh instead")
    IMPETBiCGStabILU0Solver(const Problem &problem, int overlapSize=1)
    : problem_(problem)
    , overlapSize_(overlapSize)
    {}

    ~IMPETBiCGStabILU0Solver()
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
    }

    /*!
     * \brief Actually solve the linear system of equations.
     *
     * \return true if the residual reduction could be achieved, else false.
     */
    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        int verbosityLevel = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        const int maxIter = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        const double relaxation = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, PreconditionerRelaxation);

        if (!overlapMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        }

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlapMatrix_->assignCopy(M);
        overlapb_->assignAdd(b);
        (*overlapx_) = 0.0;

        /*
        overlapMatrix_->print();
        overlapb_->print();
        exit(1);
        */

        // create sequential and overlapping preconditioners
        //SeqPreconditioner seqPreCond(*overlapMatrix_, 1, 1.0);
        SeqPreconditioner seqPreCond(*overlapMatrix_, relaxation);
        OverlappingPreconditioner preCond(seqPreCond, overlapMatrix_->overlap());

        // create the scalar products and linear operators for ISTL
        OverlappingScalarProduct scalarProd(overlapMatrix_->overlap());
        OverlappingOperator opA(*overlapMatrix_);

        // create the actual solver
        Solver solver(opA,
                      scalarProd,
                      preCond,
                      residReduction,
                      maxIter ,
                      verbosityLevel);

        // run the solver
        Dune::InverseOperatorResult result;
        solver.apply(*overlapx_, *overlapb_, result);

        // copy the result back to the non-overlapping vector
        overlapx_->assignTo(x);

        // return the result of the solver
        return result.converged;
    }

private:
    void prepare_(const Matrix &M)
    {
        ElementBorderListFromGrid<GridView, ElementMapper>
            borderListCreator(problem_.gridView(), problem_.elementMapper());

        // create the overlapping Jacobian matrix
        overlapMatrix_ = std::make_shared<OverlappingMatrix> (M,
                                                borderListCreator.foreignBorderList(),
                                                borderListCreator.domesticBorderList(),
                                                overlapSize_);

        // create the overlapping vectors for the residual and the
        // solution
        overlapb_ = std::make_shared<OverlappingVector>(overlapMatrix_->overlap());
        overlapx_ = std::make_shared<OverlappingVector>(*overlapb_);
    }

    void cleanup_()
    {
        overlapMatrix_.template reset<OverlappingMatrix>(0);
        overlapb_.template reset<OverlappingVector>(0);
        overlapx_.template reset<OverlappingVector>(0);
    }

    const Problem &problem_;

    int overlapSize_;
    std::shared_ptr<OverlappingMatrix> overlapMatrix_;
    std::shared_ptr<OverlappingVector> overlapb_;
    std::shared_ptr<OverlappingVector> overlapx_;
};

} // namespace Dumux

#endif
