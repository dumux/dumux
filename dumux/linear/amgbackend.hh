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

#if HAVE_DUNE_PDELAB

#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/backend/novlpistlsolverbackend.hh>
#include <dune/pdelab/backend/ovlpistlsolverbackend.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>

#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dumux/linear/p0fem.hh>

#include <dumux/implicit/box/boxproperties.hh>
#include <dumux/implicit/cellcentered/ccproperties.hh>
#include <dumux/decoupled/common/pressureproperties.hh>
#include "linearsolverproperties.hh"

namespace Dumux {

// forward declaration for the property definitions
template <class TypeTag> class AMGBackend;

namespace Properties
{
//! the PDELab finite element map used for the gridfunctionspace
NEW_PROP_TAG(AMGLocalFemMap);

//! box: use the (multi-)linear local FEM space associated with cubes by default
SET_PROP(BoxModel, AMGLocalFemMap)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> type;
};

//! cell-centered: use the element-wise constant local FEM space by default
SET_PROP(CCModel, AMGLocalFemMap)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef Dumux::P0LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};

//! decoupled models: use the element-wise constant local FEM space by default
SET_PROP(DecoupledModel, AMGLocalFemMap)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef Dumux::P0LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};

//! the type of the employed PDELab backend
NEW_PROP_TAG(AMGPDELabBackend);

//! box: use the non-overlapping PDELab AMG backend
SET_PROP(BoxModel, AMGPDELabBackend)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GridOperator> type;
};

//! cell-centered: use the overlapping PDELab AMG backend
SET_PROP(CCModel, AMGPDELabBackend)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GridOperator> type;
};

//! decoupled model: use the overlapping PDELab AMG backend
SET_PROP(DecoupledModel, AMGPDELabBackend)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GridOperator> type;
};

//! box: reset the type of solution vector to be PDELab conforming
SET_PROP(BoxModel, SolutionVector)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef typename GridOperator::Traits::Domain type;
};

//! cell-centered: reset the type of solution vector to be PDELab conforming
SET_PROP(CCModel, SolutionVector)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef typename GridOperator::Traits::Domain type;
};


//! decoupled model: reset the type of solution vector to be PDELab conforming
SET_PROP(DecoupledModel, PressureSolutionVector)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef typename GridOperator::Traits::Domain type;
};

//! decoupled model: reset the type of solution vector to be PDELab conforming
SET_PROP(DecoupledModel, PressureRHSVector)
{
    typedef typename Dumux::AMGBackend<TypeTag>::GridOperator GridOperator;
public:
    typedef typename GridOperator::Traits::Domain type;
};

//! set a property JacobianMatrix also for the decoupled models
SET_PROP(DecoupledModel, JacobianMatrix)
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) type;
};


}

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
    typedef typename GET_PROP_TYPE(TypeTag, AMGLocalFemMap) LocalFemMap;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { dim = GridView::dimension };
    typedef typename Dune::PDELab::NoConstraints Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef Dune::PDELab::GridFunctionSpace<GridView, 
                                            LocalFemMap, 
                                            Constraints, 
                                            Dune::PDELab::ISTLVectorBackend<numEq> 
                                           > ScalarGridFunctionSpace;
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, 
                                                 numEq, 
                                                 Dune::PDELab::GridFunctionSpaceBlockwiseMapper
                                                > GridFunctionSpace;
    typedef typename GridFunctionSpace::template ConstraintsContainer<Scalar>::Type ConstraintsTrafo;
    typedef int LocalOperator;

public:
    typedef typename Dune::PDELab::GridOperator<GridFunctionSpace,
                                       GridFunctionSpace,
                                       LocalOperator,
                                       Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
                                       Scalar, Scalar, Scalar,
                                       ConstraintsTrafo,
                                       ConstraintsTrafo,
                                       true
                                      > GridOperator;
    typedef typename GET_PROP_TYPE(TypeTag, AMGPDELabBackend) PDELabBackend;

    /*!
     * \brief Construct the backend.
     * 
     * \param problem the problem at hand
     */
    AMGBackend(const Problem& problem)
    : problem_(problem)
    {
        fem_ = Dune::make_shared<LocalFemMap>();
        constraints_ = Dune::make_shared<Constraints>();
        scalarGridFunctionSpace_ = Dune::make_shared<ScalarGridFunctionSpace>(problem.gridView(), *fem_, *constraints_);
        gridFunctionSpace_ = Dune::make_shared<GridFunctionSpace>(*scalarGridFunctionSpace_);
        int maxIt = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, MaxIterations);
        int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
        imp_ = Dune::make_shared<PDELabBackend>(*gridFunctionSpace_, maxIt, verbosity);
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
        static const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        imp_->apply(A, x, b, residReduction);

        result_.converged  = imp_->result().converged;
        result_.iterations = imp_->result().iterations;
        result_.elapsed    = imp_->result().elapsed;
        result_.reduction  = imp_->result().reduction;
        result_.conv_rate  = imp_->result().conv_rate;

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
    Dune::shared_ptr<LocalFemMap> fem_;
    Dune::shared_ptr<Constraints> constraints_;
    Dune::shared_ptr<ScalarGridFunctionSpace> scalarGridFunctionSpace_;
    Dune::shared_ptr<GridFunctionSpace> gridFunctionSpace_;
    Dune::shared_ptr<PDELabBackend> imp_;
    Dune::InverseOperatorResult result_;
};

/*!
 * \brief Provides a linear solver using the sequential PDELab AMG backend.
 */
template <class TypeTag>
class SeqAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename AMGBackend<TypeTag>::GridOperator GridOperator;
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GridOperator> PDELabBackend;
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
        imp_ = Dune::make_shared<PDELabBackend>(maxIt, verbosity);

        static const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        imp_->apply(A, x, b, residReduction);

        result_.converged  = imp_->result().converged;
        result_.iterations = imp_->result().iterations;
        result_.elapsed    = imp_->result().elapsed;
        result_.reduction  = imp_->result().reduction;
        result_.conv_rate  = imp_->result().conv_rate;

        imp_.template reset<PDELabBackend>(0);
        
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
    Dune::shared_ptr<PDELabBackend> imp_;
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
    typedef typename AMGBackend<TypeTag>::GridOperator GridOperator;
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GridOperator> PDELabBackend;
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
        imp_ = Dune::make_shared<PDELabBackend>(maxIt, verbosity);
        
        static const double residReduction = GET_PARAM_FROM_GROUP(TypeTag, double, LinearSolver, ResidualReduction);
        imp_->apply(A, x, b, residReduction);

        result_.converged  = imp_->result().converged;
        result_.iterations = imp_->result().iterations;
        result_.elapsed    = imp_->result().elapsed;
        result_.reduction  = imp_->result().reduction;
        result_.conv_rate  = imp_->result().conv_rate;

        imp_.template reset<PDELabBackend>(0);
        
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
    Dune::shared_ptr<PDELabBackend> imp_;
    Dune::InverseOperatorResult result_;
};

} // namespace Dumux

#endif // HAVE_DUNE_PDELAB
#endif // DUMUX_AMGBACKEND_HH
