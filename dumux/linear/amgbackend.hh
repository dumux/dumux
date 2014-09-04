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
 * \brief Provides linear solvers using the PDELab AMG backends.
 */
#ifndef DUMUX_AMGBACKEND_HH
#define DUMUX_AMGBACKEND_HH

#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvers.hh>

#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/linear/amgproperties.hh>
#include <dumux/linear/amgparallelhelpers.hh>
#include <dumux/linear/p0fem.hh>

namespace Dumux {

/*!
 * \brief Scale the linear system by the inverse of 
 * its (block-)diagonal entries.
 * 
 * \param matrix the matrix to scale
 * \param rhs the right hand side vector to scale
 */
template <class Matrix, class Vector>
void scaleLinearSystem(Matrix& matrix, Vector& rhs)
{
    typename Matrix::RowIterator row = matrix.begin();
    for(; row != matrix.end(); ++row)
    {
        typedef typename Matrix::size_type size_type;
        size_type rowIdx = row.index();

        typedef typename Matrix::block_type MatrixBlock;
        MatrixBlock diagonal = matrix[rowIdx][rowIdx];
        diagonal.invert();

        typedef typename Vector::block_type VectorBlock;
        const VectorBlock b = rhs[rowIdx];
        diagonal.mv(b, rhs[rowIdx]);

        typename Matrix::ColIterator col = row->begin();
        for (; col != row->end(); ++col)
            col->leftmultiply(diagonal);
    }
}

/*!
 * \brief Provides a linear solver using the parallel PDELab AMG backend.
 */
template <class TypeTag>
class AMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};

    typedef typename GET_PROP(TypeTag, AMGPDELabBackend) PDELabBackend;
    typedef typename PDELabBackend::LinearOperator LinearOperator;
    typedef typename PDELabBackend::VType VType;
    typedef typename PDELabBackend::Comm Comm;
    typedef typename PDELabBackend::Smoother Smoother;
    typedef Dune::Amg::AMG<typename PDELabBackend::LinearOperator, VType,
                           Smoother,Comm> AMGType;
    typedef typename PDELabBackend::LinearOperator::matrix_type BCRSMat;

public:    
    /*!
     * \brief Construct the backend.
     * 
     * \param problem the problem at hand
     */
    AMGBackend(const Problem& problem)
    : problem_(problem), phelper_(problem_), firstCall_(true)
    {
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
        int maxIt = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        static const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        
#if HAVE_MPI
        Dune::SolverCategory::Category category = PDELabBackend::isNonOverlapping?
            Dune::SolverCategory::nonoverlapping : Dune::SolverCategory::overlapping;

        if(PDELabBackend::isNonOverlapping && firstCall_)
        {
            phelper_.initGhostsAndOwners();
        }

        typename PDELabBackend::Comm comm(problem_.gridView().comm(), category);
        
        if(PDELabBackend::isNonOverlapping)
        {
            // extend the matrix pattern such that it is usable for AMG
            EntityExchanger<TypeTag> exchanger(problem_);
            exchanger.getExtendedMatrix(A, phelper_);
            exchanger.sumEntries(A);
        }
        phelper_.createIndexSetAndProjectForAMG(A, comm);

        typename PDELabBackend::LinearOperator fop(A, comm);
        typename PDELabBackend::ScalarProduct sp(comm);
        int rank = comm.communicator().rank();

        // Make rhs consistent
        if(PDELabBackend::isNonOverlapping)
        {
            phelper_.makeNonOverlappingConsistent(b);
        }
#else
        typename PDELabBackend::Comm  comm;
        typename PDELabBackend::LinearOperator fop(A);
        typename PDELabBackend::ScalarProduct sp;
        int rank=0;
#endif
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments
            SmootherArgs;
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<BCRSMat,
                                                                          Dune::Amg::FirstDiagonal> >
            Criterion;
        // \todo Check whether the default accumulation mode atOnceAccu is needed.
        Dune::Amg::Parameters params(15,2000,1.2,1.6,Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(GET_PROP_TYPE(TypeTag, GridView)::Traits::Grid::dimension);
        params.setDebugLevel(verbosity);
        Criterion criterion(params);
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        AMGType amg(fop, criterion, smootherArgs, comm);        
        Dune::BiCGSTABSolver<typename PDELabBackend::VType> solver(fop, sp, amg, residReduction, maxIt, 
                                                                 rank==0?verbosity: 0);


        solver.apply(x, b, result_);
        firstCall_ = false;
        return result_.converged;
        }

    /*!
     * \brief The result containing the convergence history.
     */
    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }
    
private:
    const Problem& problem_;
    ParallelISTLHelper<TypeTag> phelper_;
    Dune::InverseOperatorResult result_;
    bool firstCall_;
};

/*!
 * \brief Provides a linear solver using the sequential PDELab AMG backend.
 */
template <class TypeTag>
class SeqAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
    typedef Dune::Amg::SequentialInformation Comm;
    typedef Dune::AssembledLinearOperator<MType,VType, VType> LinearOperator;
    typedef Dune::SeqSSOR<MType,VType, VType> Smoother;
    typedef Dune::Amg::AMG<LinearOperator,VType,Smoother> AMGType;
public:

    /*!
     * \brief Construct the backend.
     * 
     * \param problem the problem at hand
     */
    SeqAMGBackend(const Problem& problem)
    : problem_(problem)
    {}

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
        int maxIt = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        static const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        LinearOperator fop(A);
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<typename LinearOperator::matrix_type,
                                                                          Dune::Amg::FirstDiagonal> >
            Criterion;
        Criterion criterion(15,2000);
        criterion.setDefaultValuesIsotropic(2);
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        AMGType amg(fop, criterion, smootherArgs, 1, 1, 1);
        Dune::BiCGSTABSolver<VType> solver(fop, amg, residReduction, maxIt, 
                                          verbosity);
        solver.apply(x, b, result_);        
        return result_.converged;
    }

    /*!
     * \brief The result containing the convergence history.
     */
    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    const Problem& problem_;
    Dune::InverseOperatorResult result_;
};

/*!
 * \brief Provides a linear solver using the sequential PDELab AMG backend.
 * 
 * The linear system is scaled beforehand, possibly improving the 
 * convergence behavior of the iterative solver.
 */
template <class TypeTag>
class ScaledSeqAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > MType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,numEq> > VType;
    typedef Dune::Amg::SequentialInformation Comm;
    typedef Dune::MatrixAdapter<MType,VType, VType> LinearOperator;
    typedef Dune::SeqSSOR<MType,VType, VType> Smoother;
    typedef Dune::Amg::AMG<LinearOperator,VType,Smoother> AMGType;
public:

    /*!
     * \brief Construct the backend.
     * 
     * \param problem the problem at hand
     */
    ScaledSeqAMGBackend(const Problem& problem)
    : problem_(problem)
    {}

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
        scaleLinearSystem(A, b);

        int maxIt = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        static const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        LinearOperator fop(A);
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MType,
                                                                          Dune::Amg::FirstDiagonal> >
            Criterion;
        Criterion criterion(15,2000);
        criterion.setDefaultValuesIsotropic(2);
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;
        AMGType amg(fop, criterion, smootherArgs, 1, 1, 1);
        Dune::BiCGSTABSolver<VType> solver(fop, amg, residReduction, maxIt, 
                                          verbosity);
        
        return result_.converged;
    }

    /*!
     * \brief The result containing the convergence history.
     */
    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    const Problem& problem_;
    Dune::InverseOperatorResult result_;
};

} // namespace Dumux

#endif // DUMUX_AMGBACKEND_HH
