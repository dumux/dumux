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
 * \ingroup Linear
 *
 * \brief Provides a linear solver based on the ISTL AMG.
 */
#ifndef DUMUX_AMGBACKEND_HH
#define DUMUX_AMGBACKEND_HH

#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvers.hh>

#include <dumux/linear/solver.hh>
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/linear/amgproperties.hh>
#include <dumux/linear/amgparallelhelpers.hh>

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
 * \brief Provides a linear solver using the ISTL AMG.
 */
template <class TypeTag>
class AMGBackend : public LinearSolver<TypeTag>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using AmgTraits = typename GET_PROP(TypeTag, AmgTraits);
    enum { numEq = AmgTraits::numEq };
    using LinearOperator = typename AmgTraits::LinearOperator;
    using ScalarProduct = typename AmgTraits::ScalarProduct;
    using VType = typename AmgTraits::VType;
    using Comm = typename AmgTraits::Comm;
    using Smoother = typename AmgTraits::Smoother;
    using AMGType = Dune::Amg::AMG<typename AmgTraits::LinearOperator, VType, Smoother,Comm>;
    using BCRSMat = typename AmgTraits::LinearOperator::matrix_type;
    using DofMapper = typename AmgTraits::DofMapper;

    const std::string modelParamGroup = GET_PROP_VALUE(TypeTag, ModelParameterGroup);

public:
    /*!
     * \brief Construct the backend.
     *
     * \param problem the problem at hand
     */
    AMGBackend(const GridView& gridView, const DofMapper& mapper)
    : LinearSolver<TypeTag>(gridView, mapper)
    , phelper_(gridView, mapper)
    , firstCall_(true)
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
        const int verbosity = getParamFromGroup<int>(modelParamGroup, "LinearSolver.Verbosity");
        const int maxIter = getParamFromGroup<int>(modelParamGroup, "LinearSolver.MaxIterations");
        const double residReduction = getParamFromGroup<double>(modelParamGroup, "LinearSolver.ResidualReduction");

        int rank = 0;
        std::shared_ptr<Comm> comm;
        std::shared_ptr<LinearOperator> fop;
        std::shared_ptr<ScalarProduct> sp;
        static const int dofCodim = AmgTraits::dofCodim;
        static const bool isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v;
        prepareLinearAlgebra_<Matrix, Vector, isParallel>(A, b, rank, comm, fop, sp);

        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<BCRSMat, Dune::Amg::FirstDiagonal>>;

        //! \todo Check whether the default accumulation mode atOnceAccu is needed.
        Dune::Amg::Parameters params(15,2000,1.2,1.6,Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(Grid::dimension);
        params.setDebugLevel(verbosity);
        Criterion criterion(params);
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        AMGType amg(*fop, criterion, smootherArgs, *comm);
        Dune::BiCGSTABSolver<VType> solver(*fop, *sp, amg, residReduction, maxIter,
                                           rank == 0 ? verbosity : 0);

        solver.apply(x, b, result_);
        firstCall_ = false;
        return result_.converged;
    }

    std::string name() const
    {
        return "AMG preconditioned BiCGSTAB solver";
    }

    /*!
     * \brief The result containing the convergence history.
     */
    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:

    /*!
     * \brief Prepare the linear algebra member variables.
     *
     * At compile time, correct constructor calls have to be chosen,
     * depending on whether the setting is parallel or sequential.
     * Since several template parameters are present, this cannot be solved
     * by a full function template specialization. Instead, class template
     * specialization has to be used.
     * This adapts example 4 from http://www.gotw.ca/publications/mill17.htm.

     * The function is called from the solve function. The call is
     * forwarded to the corresponding function of a class template.
     *
     * \tparam isParallel decides if the setting is parallel or sequential
     */
    template<class Matrix, class Vector, bool isParallel>
    void prepareLinearAlgebra_(Matrix& A, Vector& b, int& rank,
                               std::shared_ptr<Comm>& comm,
                               std::shared_ptr<LinearOperator>& fop,
                               std::shared_ptr<ScalarProduct>& sp)
    {
        LinearAlgebraPreparator<TypeTag, isParallel>
          ::prepareLinearAlgebra(A, b, rank, comm, fop, sp,
                                 this->gridView(), this->dofMapper(), phelper_, firstCall_);
    }

    ParallelISTLHelper<TypeTag> phelper_;
    Dune::InverseOperatorResult result_;
    bool firstCall_;
};

} // namespace Dumux

#endif // DUMUX_AMGBACKEND_HH
