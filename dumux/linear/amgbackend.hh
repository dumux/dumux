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
 * \brief Provides a linear solver using the PDELab AMG backend.
 */
#ifndef DUMUX_AMGBACKEND_HH
#define DUMUX_AMGBACKEND_HH

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include "linearsolverproperties.hh"
#include "novlpistlsolverbackend.hh"
#include "ovlpistlsolverbackend.hh"
#include "seqistlsolverbackend.hh"

namespace Dumux {

namespace Properties
{
NEW_PROP_TAG(GridOperator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ImplicitIsBox);
NEW_PROP_TAG(NumEq);
}

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
 * \brief Provides a linear solver using the PDELab AMG backend.
 */
template <class TypeTag, class LocalFemMap, class PDELabBackend>
class AMGBackend
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { dim = GridView::dimension };

    typedef typename Dune::PDELab::NoConstraints Constraints;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::PDELab::GridFunctionSpace<GridView, 
                                            LocalFemMap, 
                                            Constraints, 
                                            Dumux::ISTLVectorBackend<numEq> 
                                           > ScalarGridFunctionSpace;
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, 
                                                 numEq, 
                                                 Dune::PDELab::GridFunctionSpaceBlockwiseMapper
                                                > GridFunctionSpace;
    typedef typename GridFunctionSpace::template ConstraintsContainer<Scalar>::Type ConstraintsTrafo;
    typedef int LocalOperator;

public:
    typedef Dune::PDELab::GridOperator<GridFunctionSpace,
                                       GridFunctionSpace,
                                       LocalOperator,
                                       Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
                                       Scalar, Scalar, Scalar,
                                       ConstraintsTrafo,
                                       ConstraintsTrafo,
                                       true
                                      > GridOperator;

    AMGBackend(const Problem& problem)
    : problem_(problem)
    {}

    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        imp_ = new PDELabBackend(problem_.gridFunctionSpace(),
                GET_PROP_VALUE(TypeTag, LinearSolverMaxIterations),
                GET_PROP_VALUE(TypeTag, LinearSolverVerbosity));

        static const double residReduction = GET_PROP_VALUE(TypeTag, LinearSolverResidualReduction);
        imp_->apply(A, x, b, residReduction);

        result_.converged  = imp_->result().converged;
        result_.iterations = imp_->result().iterations;
        result_.elapsed    = imp_->result().elapsed;
        result_.reduction  = imp_->result().reduction;
        result_.conv_rate  = imp_->result().conv_rate;

        delete imp_;

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    const Problem& problem_;
    PDELabBackend *imp_;
    Dune::InverseOperatorResult result_;
};

template <class TypeTag>
class CCAMGBackend : public AMGBackend<TypeTag, 
                                       Dune::PDELab::P0LocalFiniteElementMap<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                             typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                             AMGBackend::GridView::dimension>, 
                                       ISTLBackend_BCGS_AMG_SSOR<typename AMGBackend::GridOperator> >
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename AMGBackend<TypeTag, 
                                Dune::PDELab::P0LocalFiniteElementMap<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                      typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                      typename AMGBackend::GridView::dimension>, 
                                ISTLBackend_BCGS_AMG_SSOR<typename AMGBackend::GridOperator> > ParentType;
public:
    CCAMGBackend(const Problem& problem)
    : ParentType(problem)
    {}
};

template <class TypeTag>
class BoxAMGBackend : public AMGBackend<TypeTag, 
                                       Dune::PDELab::Q1LocalFiniteElementMap<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                             typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                             typename AMGBackend::GridView::dimension>, 
                                       ISTLBackend_NOVLP_BCGS_AMG_SSOR<typename AMGBackend::GridOperator> >
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename AMGBackend<TypeTag, 
                                Dune::PDELab::P0LocalFiniteElementMap<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                      typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                      typename AMGBackend::GridView::dimension>, 
                                ISTLBackend_NOVLP_BCGS_AMG_SSOR<typename AMGBackend::GridOperator> > ParentType;
public:
    BoxAMGBackend(const Problem& problem)
    : ParentType(problem)
    {}
};

template <class TypeTag>
class SeqAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef ISTLBackend_SEQ_BCGS_AMG_SSOR<GridOperator> PDELabBackend;
public:

    SeqAMGBackend(const Problem& problem)
    : problem_(problem)
    {}

    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        imp_ = new PDELabBackend(
                GET_PROP_VALUE(TypeTag, LinearSolverMaxIterations),
                GET_PROP_VALUE(TypeTag, LinearSolverVerbosity));

        static const double residReduction = GET_PROP_VALUE(TypeTag, LinearSolverResidualReduction);
        imp_->apply(A, x, b, residReduction);

        result_.converged  = imp_->result().converged;
        result_.iterations = imp_->result().iterations;
        result_.elapsed    = imp_->result().elapsed;
        result_.reduction  = imp_->result().reduction;
        result_.conv_rate  = imp_->result().conv_rate;

        delete imp_;

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    const Problem& problem_;
    PDELabBackend *imp_;
    Dune::InverseOperatorResult result_;
};

template <class TypeTag>
class ScaledSeqAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef ISTLBackend_SEQ_BCGS_AMG_SSOR<GridOperator> PDELabBackend;
public:

    ScaledSeqAMGBackend(const Problem& problem)
    : problem_(problem)
    {}

    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        scaleLinearSystem(A, b);

        imp_ = new PDELabBackend(
                GET_PROP_VALUE(TypeTag, LinearSolverMaxIterations),
                GET_PROP_VALUE(TypeTag, LinearSolverVerbosity));

        static const double residReduction = GET_PROP_VALUE(TypeTag, LinearSolverResidualReduction);
        imp_->apply(A, x, b, residReduction);

        result_.converged  = imp_->result().converged;
        result_.iterations = imp_->result().iterations;
        result_.elapsed    = imp_->result().elapsed;
        result_.reduction  = imp_->result().reduction;
        result_.conv_rate  = imp_->result().conv_rate;

        delete imp_;

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    const Problem& problem_;
    PDELabBackend *imp_;
    Dune::InverseOperatorResult result_;
};

} // namespace Dumux

#endif
