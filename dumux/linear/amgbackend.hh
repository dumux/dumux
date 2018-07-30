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
 * \brief Provides a parallel linear solver based on the ISTL AMG preconditioner
 *        and the ISTL BiCGSTAB solver.
 */
#ifndef DUMUX_PARALLEL_AMGBACKEND_HH
#define DUMUX_PARALLEL_AMGBACKEND_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/solvers.hh>

#include <dumux/linear/solver.hh>
#include <dumux/linear/amgparallelhelpers.hh>

namespace Dumux {

/*!
 * \ingroup Linear
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
        using size_type = typename Matrix::size_type;
        size_type rowIdx = row.index();

        using MatrixBlock = typename Matrix::block_type;
        MatrixBlock diagonal = matrix[rowIdx][rowIdx];
        diagonal.invert();

        using VectorBlock = typename Vector::block_type;
        const VectorBlock b = rhs[rowIdx];
        diagonal.mv(b, rhs[rowIdx]);

        typename Matrix::ColIterator col = row->begin();
        for (; col != row->end(); ++col)
            col->leftmultiply(diagonal);
    }
}

/*!
 * \ingroup Linear
 * \brief A linear solver based on the ISTL AMG preconditioner
 *        and the ISTL BiCGSTAB solver.
 */
template <class GridView, class AmgTraits>
class ParallelAMGBackend : public LinearSolver
{
    using Grid = typename GridView::Grid;
    using LinearOperator = typename AmgTraits::LinearOperator;
    using ScalarProduct = typename AmgTraits::ScalarProduct;
    using VType = typename AmgTraits::VType;
    using Comm = typename AmgTraits::Comm;
    using Smoother = typename AmgTraits::Smoother;
    using AMGType = Dune::Amg::AMG<typename AmgTraits::LinearOperator, VType, Smoother,Comm>;
    using BCRSMat = typename AmgTraits::LinearOperator::matrix_type;
    using DofMapper = typename AmgTraits::DofMapper;
public:
    /*!
     * \brief Construct the backend for the sequential case only
     *
     * \param gridView the grid view on which we are performing the multi-grid
     * \param dofMapper an index mapper for dof entities
     * \param paramGroup the parameter group for parameter lookup
     */
    ParallelAMGBackend(const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , firstCall_(true)
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
            DUNE_THROW(Dune::InvalidStateException, "Using sequential constructor for parallel run. Use signature with gridView and dofMapper!");
    }

    /*!
     * \brief Construct the backend for parallel or sequential runs
     *
     * \param gridView the grid view on which we are performing the multi-grid
     * \param dofMapper an index mapper for dof entities
     * \param paramGroup the parameter group for parameter lookup
     */
    ParallelAMGBackend(const GridView& gridView,
                       const DofMapper& dofMapper,
                       const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , phelper_(std::make_shared<ParallelISTLHelper<GridView, AmgTraits>>(gridView, dofMapper))
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
        int rank = 0;
        std::shared_ptr<Comm> comm;
        std::shared_ptr<LinearOperator> fop;
        std::shared_ptr<ScalarProduct> sp;
        static const int dofCodim = AmgTraits::dofCodim;
        static const bool isParallel = AmgTraits::isParallel;
        prepareLinearAlgebra_<Matrix, Vector, isParallel>(A, b, rank, comm, fop, sp);

        using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
        using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<BCRSMat, Dune::Amg::FirstDiagonal>>;

        //! \todo Check whether the default accumulation mode atOnceAccu is needed.
        //! \todo make parameters changeable at runtime from input file / parameter tree
        Dune::Amg::Parameters params(15,2000,1.2,1.6,Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(Grid::dimension);
        params.setDebugLevel(this->verbosity());
        Criterion criterion(params);
        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1;

        AMGType amg(*fop, criterion, smootherArgs, *comm);
        Dune::BiCGSTABSolver<VType> solver(*fop, *sp, amg, this->residReduction(), this->maxIter(),
                                           rank == 0 ? this->verbosity() : 0);

        solver.apply(x, b, result_);
        firstCall_ = false;
        return result_.converged;
    }

    /*!
     * \brief The name of the solver
     */
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
     * \tparam Matrix the matrix type
     * \tparam Vector the vector type
     * \tparam isParallel decides if the setting is parallel or sequential
     */
    template<class Matrix, class Vector, bool isParallel>
    void prepareLinearAlgebra_(Matrix& A, Vector& b, int& rank,
                               std::shared_ptr<Comm>& comm,
                               std::shared_ptr<LinearOperator>& fop,
                               std::shared_ptr<ScalarProduct>& sp)
    {
        LinearAlgebraPreparator<GridView, AmgTraits, isParallel>
          ::prepareLinearAlgebra(A, b, rank, comm, fop, sp,
                                 *phelper_, firstCall_);
    }

    std::shared_ptr<ParallelISTLHelper<GridView, AmgTraits>> phelper_;
    Dune::InverseOperatorResult result_;
    bool firstCall_;
};

} // namespace Dumux

#include <dumux/common/properties.hh>
#include <dumux/linear/amgtraits.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A linear solver based on the ISTL AMG preconditioner
 *        and the ISTL BiCGSTAB solver.
 * \note This is an adaptor using a TypeTag
 */
template<class TypeTag>
using AMGBackend = ParallelAMGBackend<typename GET_PROP_TYPE(TypeTag, GridView), AmgTraits<typename GET_PROP_TYPE(TypeTag, JacobianMatrix),
                                                                                           typename GET_PROP_TYPE(TypeTag, SolutionVector),
                                                                                           typename GET_PROP_TYPE(TypeTag, FVGridGeometry)>>;

} // namespace Dumux

#endif // DUMUX_AMGBACKEND_HH
