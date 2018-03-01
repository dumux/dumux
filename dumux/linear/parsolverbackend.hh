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
 * \brief Dumux parallel linear solver backends
 */
#ifndef DUMUX_PAR_SOLVER_BACKEND_HH
#define DUMUX_PAR_SOLVER_BACKEND_HH

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dumux/linear/solver.hh>

// the AMG-BiCGSTABBackend
#include <dumux/linear/amgbackend.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A linear solver based on the ISTL SSOR preconditioner
 *        and the ISTL CG solver.
 */
template <class TypeTag>
class ParallelSSORCGBackend : public LinearSolver
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using Traits = AmgTraits<TypeTag>;
    using Grid = typename GridView::Grid;
    enum { numEq = Traits::numEq };
    using LinearOperator = typename Traits::LinearOperator;
    using ScalarProduct = typename Traits::ScalarProduct;
    using VType = typename Traits::VType;
    using Comm = typename Traits::Comm;
    using BCRSMat = typename Traits::LinearOperator::matrix_type;
    using DofMapper = typename Traits::DofMapper;
    using SeqPreconditioner = Dune::SeqSSOR<BCRSMat, VType, VType>;
    using Preconditioner = Dune::BlockPreconditioner<VType, VType, Comm, SeqPreconditioner>;
public:
    /*!
     * \brief Construct the backend for the sequential case only
     *
     * \param gridView the grid view on which we are performing the multi-grid
     * \param dofMapper an index mapper for dof entities
     * \param paramGroup the parameter group for parameter lookup
     */
    ParallelSSORCGBackend(const std::string& paramGroup = "")
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
    ParallelSSORCGBackend(const GridView& gridView,
                          const DofMapper& dofMapper,
                          const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , phelper_(std::make_shared<ParallelISTLHelper<GridView, Traits>>(gridView, dofMapper))
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
        static const int dofCodim = Traits::dofCodim;
        static const bool isParallel = Dune::Capabilities::canCommunicate<Grid, dofCodim>::v;
        prepareLinearAlgebra_<Matrix, Vector, isParallel>(A, b, rank, comm, fop, sp);

        SeqPreconditioner seqSSOR(A, this->precondIter(), this->relaxation());
        Preconditioner parSSOR(seqSSOR, *comm);
        Dune::CGSolver<VType> solver(*fop, *sp, parSSOR, this->residReduction(), this->maxIter(),
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
        return "SSOR preconditioned CG solver";
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
        LinearAlgebraPreparator<GridView, Traits, isParallel>
          ::prepareLinearAlgebra(A, b, rank, comm, fop, sp,
                                 *phelper_, firstCall_);
    }

    std::shared_ptr<ParallelISTLHelper<GridView, Traits>> phelper_;
    Dune::InverseOperatorResult result_;
    bool firstCall_;
};

} // end namespace Dumux

#endif
