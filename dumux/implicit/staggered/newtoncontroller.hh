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
 * \brief A 2p1cni specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_STAGGERED_NEWTON_CONTROLLER_HH
#define DUMUX_STAGGERED_NEWTON_CONTROLLER_HH

#include "properties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#include "newtonconvergencewriter.hh"

namespace Dumux {

namespace Properties
{
    SET_PROP(StaggeredModel, LinearSolverBlockSize)
    {
        // LinearSolverAcceptsMultiTypeMatrix<T>::value
        // TODO: make somehow dependend? or only relevant for direct solvers?
    public:
        static constexpr auto value = 1;
    };
}
/*!
 * \ingroup PNMModel
 * \brief A PNM specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class StaggeredNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef NewtonConvergenceWriter<TypeTag> StaggeredNewtonConvergenceWriter;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
        numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter),
        numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace)
    };

public:
    StaggeredNewtonController(const Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * If the linear solver doesn't accept multitype matrices we copy the matrix
     * into a 1x1 block BCRS matrix for solving.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<typename T = TypeTag>
    typename std::enable_if<!LinearSolverAcceptsMultiTypeMatrix<T>::value, void>::type
    newtonSolveLinear(JacobianMatrix &A,
                      SolutionVector &x,
                      SolutionVector &b)
    {
        try
        {
            if (this->numSteps_ == 0)
                this->initialResidual_ = b.two_norm();

            // copy the matrix and the vector to types the IterativeSolverBackend can handle
            using MatrixBlock = typename Dune::FieldMatrix<Scalar, 1, 1>;
            using SparseMatrix = typename Dune::BCRSMatrix<MatrixBlock>;

            // get the new matrix sizes
            std::size_t numRows = numEqCellCenter*A[cellCenterIdx][cellCenterIdx].N() + numEqFace*A[faceIdx][cellCenterIdx].N();
            std::size_t numCols = numEqCellCenter*A[cellCenterIdx][cellCenterIdx].M() + numEqFace*A[cellCenterIdx][faceIdx].M();

            // check matrix sizes
            assert(A[cellCenterIdx][cellCenterIdx].N() == A[cellCenterIdx][faceIdx].N());
            assert(A[faceIdx][cellCenterIdx].N() == A[faceIdx][faceIdx].N());
            assert(numRows == numCols);

            // create the bcrs matrix the IterativeSolver backend can handle
            auto M = SparseMatrix(numRows, numCols, SparseMatrix::random);

            // set the rowsizes
            // A11 and A12
            for (auto row = A[cellCenterIdx][cellCenterIdx].begin(); row != A[cellCenterIdx][cellCenterIdx].end(); ++row)
                for (std::size_t i = 0; i < numEqCellCenter; ++i)
                    M.setrowsize(numEqCellCenter*row.index() + i, row->size()*numEqCellCenter);
            for (auto row = A[cellCenterIdx][faceIdx].begin(); row != A[cellCenterIdx][faceIdx].end(); ++row)
                for (std::size_t i = 0; i < numEqCellCenter; ++i)
                    M.setrowsize(numEqCellCenter*row.index() + i, M.getrowsize(numEqCellCenter*row.index() + i) + row->size()*numEqFace);
            // A21 and A22
            for (auto row = A[faceIdx][cellCenterIdx].begin(); row != A[faceIdx][cellCenterIdx].end(); ++row)
                for (std::size_t i = 0; i < numEqFace; ++i)
                    M.setrowsize(numEqFace*row.index() + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, row->size()*numEqCellCenter);
            for (auto row = A[faceIdx][faceIdx].begin(); row != A[faceIdx][faceIdx].end(); ++row)
                for (std::size_t i = 0; i < numEqFace; ++i)
                    M.setrowsize(numEqFace*row.index() + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, M.getrowsize(numEqFace*row.index() + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter) + row->size()*numEqFace);
            M.endrowsizes();

            // set the indices
            for (auto row = A[cellCenterIdx][cellCenterIdx].begin(); row != A[cellCenterIdx][cellCenterIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqCellCenter; ++i)
                        for (std::size_t j = 0; j < numEqCellCenter; ++j)
                            M.addindex(row.index()*numEqCellCenter + i, col.index()*numEqCellCenter + j);

            for (auto row = A[cellCenterIdx][faceIdx].begin(); row != A[cellCenterIdx][faceIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqCellCenter; ++i)
                        for (std::size_t j = 0; j < numEqFace; ++j)
                            M.addindex(row.index()*numEqCellCenter + i, col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter);

            for (auto row = A[faceIdx][cellCenterIdx].begin(); row != A[faceIdx][cellCenterIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqFace; ++i)
                        for (std::size_t j = 0; j < numEqCellCenter; ++j)
                            M.addindex(row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, col.index()*numEqCellCenter + j);

            for (auto row = A[faceIdx][faceIdx].begin(); row != A[faceIdx][faceIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqFace; ++i)
                        for (std::size_t j = 0; j < numEqFace; ++j)
                            M.addindex(row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter);
            M.endindices();

            // copy values
            for (auto row = A[cellCenterIdx][cellCenterIdx].begin(); row != A[cellCenterIdx][cellCenterIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqCellCenter; ++i)
                        for (std::size_t j = 0; j < numEqCellCenter; ++j)
                            M[row.index()*numEqCellCenter + i][col.index()*numEqCellCenter + j] = A[cellCenterIdx][cellCenterIdx][row.index()][col.index()][i][j];

            for (auto row = A[cellCenterIdx][faceIdx].begin(); row != A[cellCenterIdx][faceIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqCellCenter; ++i)
                        for (std::size_t j = 0; j < numEqFace; ++j)
                            M[row.index()*numEqCellCenter + i][col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter] = A[cellCenterIdx][faceIdx][row.index()][col.index()][i][j];

            for (auto row = A[faceIdx][cellCenterIdx].begin(); row != A[faceIdx][cellCenterIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqFace; ++i)
                        for (std::size_t j = 0; j < numEqCellCenter; ++j)
                            M[row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter][col.index()*numEqCellCenter + j] = A[faceIdx][cellCenterIdx][row.index()][col.index()][i][j];

            for (auto row = A[faceIdx][faceIdx].begin(); row != A[faceIdx][faceIdx].end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < numEqFace; ++i)
                        for (std::size_t j = 0; j < numEqFace; ++j)
                            M[row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter][col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter] = A[faceIdx][faceIdx][row.index()][col.index()][i][j];

            // create the vector the IterativeSolver backend can handle
            using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
            using BlockVector = typename Dune::BlockVector<VectorBlock>;

            BlockVector y, bTmp;
            y.resize(numRows);
            bTmp.resize(numCols);
            for (std::size_t i = 0; i < b[cellCenterIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    bTmp[i*numEqCellCenter + j] = b[cellCenterIdx][i][j];
            for (std::size_t i = 0; i < b[faceIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    bTmp[i*numEqFace + j + b[cellCenterIdx].N()*numEqCellCenter] = b[faceIdx][i][j];

            // printmatrix(std::cout, M, "", "");

            // testing stuff:
            BlockVector yTestUMFPack;
            yTestUMFPack.resize(A[faceIdx][faceIdx].N());
            auto yTestBICGSTAB = yTestUMFPack;
            auto yTestgmres = yTestUMFPack;
            auto yTestSUPERLU = yTestUMFPack;
            auto yTestILUPack = yTestUMFPack;
            const auto bTmpTest = b[faceIdx];

            Dune::Timer timer;
            std::cout << std::endl;

            UMFPackBackend<TypeTag> umfSolver(this->problem_());
            umfSolver.solve(A[faceIdx][faceIdx], yTestUMFPack, bTmpTest);
            std::cout << "solving velo with UMFPackBackend took " << timer.elapsed() << std::endl;
            timer.reset();

            auto getDiff = [&yTestUMFPack] (const auto& v)
            {
                Scalar diff = 0.0;

                for(int i = 0; i < yTestUMFPack.size(); ++i)
                    diff += (yTestUMFPack[i] - v[i])*(yTestUMFPack[i] - v[i]);

                return diff;
            };

            ILUPackBackend<TypeTag> iluSolver(this->problem_());
            iluSolver.solve(A[faceIdx][faceIdx], yTestILUPack, bTmpTest);
            std::cout << "solving velo with ILUPackBackend took " << timer.elapsed() << "; difference to UMFPACK: "  << getDiff(yTestILUPack) << std::endl;
            timer.reset();

            ILU0BiCGSTABBackend<TypeTag> bicgstabSolver(this->problem_());
            bicgstabSolver.solve(A[faceIdx][faceIdx], yTestBICGSTAB, bTmpTest);
            std::cout << "solving velo with ILU0BiCGSTABBackend took " << timer.elapsed()  << "; difference to UMFPACK: "  << getDiff(yTestBICGSTAB) << std::endl;
            timer.reset();

            ILU0RestartedGMResBackend<TypeTag> gmresSolver(this->problem_());
            gmresSolver.solve(A[faceIdx][faceIdx], yTestgmres, bTmpTest);
            std::cout << "solving velo with ILU0RestartedGMResBackend took " << timer.elapsed() << "; difference to UMFPACK: "  << getDiff(yTestgmres) << std::endl;
            timer.reset();


            SuperLUBackend<TypeTag> superluSolver(this->problem_());
            superluSolver.solve(A[faceIdx][faceIdx], yTestSUPERLU, bTmpTest);
            std::cout << "solving velo with SuperLUBackend took " << timer.elapsed() << "; difference to UMFPACK: "  << getDiff(yTestSUPERLU) << std::endl;
            timer.reset();


            // solve
            bool converged = this->linearSolver_.solve(M, y, bTmp);

            // copy back the result y into x
            for (std::size_t i = 0; i < x[cellCenterIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    x[cellCenterIdx][i][j] = y[i*numEqCellCenter + j];
            for (std::size_t i = 0; i < x[faceIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    x[faceIdx][i][j] = y[i*numEqFace + j + x[cellCenterIdx].N()*numEqCellCenter];

            if (!converged)
                DUNE_THROW(NumericalProblem, "Linear solver did not converge");
        }
        catch (const Dune::Exception &e)
        {
            Dumux::NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<typename T = TypeTag>
    typename std::enable_if<LinearSolverAcceptsMultiTypeMatrix<T>::value, void>::type
    newtonSolveLinear(JacobianMatrix &A,
                      SolutionVector &x,
                      SolutionVector &b)
    {
        try
        {
            if (this->numSteps_ == 0)
                this->initialResidual_ = b.two_norm();

            bool converged = this->linearSolver_.solve(A, x, b);

            if (!converged)
                DUNE_THROW(NumericalProblem, "Linear solver did not converge");
        }
        catch (const Dune::Exception &e)
        {
            Dumux::NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

     /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        if (this->enableShiftCriterion_)
            this->newtonUpdateShift(uLastIter, deltaU);

        this->writeConvergence_(uLastIter, deltaU);

        if (this->useLineSearch_)
        {
            this->lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        }
        else {
            for (unsigned int i = 0; i < uLastIter[cellCenterIdx].size(); ++i) {
                uCurrentIter[cellCenterIdx][i] = uLastIter[cellCenterIdx][i];
                uCurrentIter[cellCenterIdx][i] -= deltaU[cellCenterIdx][i];
            }
            for (unsigned int i = 0; i < uLastIter[faceIdx].size(); ++i) {
                uCurrentIter[faceIdx][i] = uLastIter[faceIdx][i];
                uCurrentIter[faceIdx][i] -= deltaU[faceIdx][i];
            }

            if (this->enableResidualCriterion_)
            {
                SolutionVector tmp(uLastIter);
                this->reduction_ = this->method().model().globalResidual(tmp, uCurrentIter);
                this->reduction_ /= this->initialResidual_;
            }
        }
    }

     /*!
     * \brief Update the maximum relative shift of the solution compared to
     *        the previous iteration.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    void newtonUpdateShift(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        this->shift_ = 0;

        for (int i = 0; i < int(uLastIter[cellCenterIdx].size()); ++i) {
            auto uNewI = uLastIter[cellCenterIdx][i];
            uNewI -= deltaU[cellCenterIdx][i];

            Scalar shiftAtDof = this->model_().relativeShiftAtDof(uLastIter[cellCenterIdx][i],
                                                            uNewI);
            this->shift_ = std::max(this->shift_, shiftAtDof);
        }
        for (int i = 0; i < int(uLastIter[faceIdx].size()); ++i) {
            auto uNewI = uLastIter[faceIdx][i];
            uNewI -= deltaU[faceIdx][i];

            Scalar shiftAtDof = this->model_().relativeShiftAtDof(uLastIter[faceIdx][i],
                                                            uNewI);
            this->shift_ = std::max(this->shift_, shiftAtDof);
        }

        if (this->gridView_().comm().size() > 1)
            this->shift_ = this->gridView_().comm().max(this->shift_);
    }

};
}

#endif
