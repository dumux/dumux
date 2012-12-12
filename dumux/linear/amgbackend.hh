/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_AMGBACKEND_HH
#define DUMUX_AMGBACKEND_HH

#include <dumux/linear/boxlinearsolver.hh>
#include <dumux/boxmodels/pdelab/pdelabadapter.hh>
#include <dumux/boxmodels/pdelab/pdelabboxistlvectorbackend.hh>
#include <dune/pdelab/backend/novlpistlsolverbackend.hh>
#include <dune/pdelab/backend/ovlpistlsolverbackend.hh>
#include <dune/pdelab/backend/seqistlsolverbackend.hh>

namespace Dumux {

namespace Properties
{
NEW_PROP_TAG(GridOperator);
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

template <class TypeTag>
class AMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_AMG_SSOR<GridOperator> Imp;
public:

    AMGBackend(const Problem& problem)
    : problem_(problem)
    {}

    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        imp_ = new Imp(problem_.gridFunctionSpace(),
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
    Imp *imp_;
    Dune::InverseOperatorResult result_;
};

template <class TypeTag>
class CCAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GridOperator> Imp;
public:

    CCAMGBackend(const Problem& problem)
    : problem_(problem)
    {}

    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        imp_ = new Imp(problem_.gridFunctionSpace(),
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
    Imp *imp_;
    Dune::InverseOperatorResult result_;
};

template <class TypeTag>
class SeqAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GridOperator> Imp;
public:

    SeqAMGBackend(const Problem& problem)
    : problem_(problem)
    {}

    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        imp_ = new Imp(
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
    Imp *imp_;
    Dune::InverseOperatorResult result_;
};

template <class TypeTag>
class ScaledSeqAMGBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GridOperator> Imp;
public:

    ScaledSeqAMGBackend(const Problem& problem)
    : problem_(problem)
    {}

    template<class Matrix, class Vector>
    bool solve(Matrix& A, Vector& x, Vector& b)
    {
        scaleLinearSystem(A, b);

        imp_ = new Imp(
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
    Imp *imp_;
    Dune::InverseOperatorResult result_;
};

} // namespace Dumux

#endif
