// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Krylov solvers with the p-multigrid preconditioner for higher-order Lagrange discretizations
 */
#ifndef DUMUX_LINEAR_PMULTIGRID_SOLVERS_HH
#define DUMUX_LINEAR_PMULTIGRID_SOLVERS_HH

#include <memory>
#include <string>
#include <utility>
#include <variant>

#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/istl/solvers.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/pmultigridpreconditioner.hh>

namespace Dumux::Detail {

/*!
 * \ingroup Linear
 * \brief An iterative dune-istl solver with the p-multigrid preconditioner (PMultigridPreconditionerFactory)
 *
 * Constructed from the grid geometry of a higher-order Lagrange discretization, which defines the prolongation
 * from the P1/Q1 space on the grid vertices. Runs sequentially and in parallel on nonoverlapping and overlapping
 * decompositions.
 */
template<class LinearSolverTraits, class LinearAlgebraTraits, class InverseOperator>
class PMultigridIstlSolver
: public IstlIterativeLinearSolver<LinearSolverTraits, LinearAlgebraTraits, InverseOperator,
                                   PMultigridPreconditionerFactory<typename LinearSolverTraits::GridView,
                                                                   typename LinearAlgebraTraits::Matrix>>
{
    using Matrix = typename LinearAlgebraTraits::Matrix;
    using Vector = typename LinearAlgebraTraits::Vector;
    using Factory = PMultigridPreconditionerFactory<typename LinearSolverTraits::GridView, Matrix>;
    using ParentType = IstlIterativeLinearSolver<LinearSolverTraits, LinearAlgebraTraits, InverseOperator, Factory>;
    static_assert(!isMultiTypeBlockMatrix<Matrix>::value, "The p-multigrid preconditioner needs a Dune::BCRSMatrix");

public:
    using ParameterInitializer = std::variant<std::string, Dune::ParameterTree>;

    /*!
     * \brief Constructor
     * \param gridGeometry the grid geometry of the discretization (FE or hybrid CVFE PQk)
     * \param params a parameter group (parameters read from the input file) or a parameter tree
     */
    template<class GridGeometry>
    explicit PMultigridIstlSolver(const GridGeometry& gridGeometry, const ParameterInitializer& params = "")
    : ParentType(gridGeometry.gridView(), gridGeometry.dofMapper(), params)
    { this->setPreconditionerFactory(Factory(gridGeometry)); }

    using ParentType::solve;

    /*!
     * \brief Solve the linear system Ax = b
     * \note A is the matrix as assembled; in a parallel run, the solver prepares it (sums its border rows).
     *       A matrix set with setMatrix() stays the one the solver reuses.
     */
    IstlSolverResult solve(Matrix& A, Vector& x, Vector& b)
    {
        // the stored solver is reconstructed (e.g. by setResidualReduction()) with the factory of the stored matrix
        struct RestoreFactory
        {
            Factory& factory;
            Factory factoryOfStoredMatrix;
            ~RestoreFactory() { factory = std::move(factoryOfStoredMatrix); }
        } restore{this->preconditionerFactory(), this->preconditionerFactory()};

        this->preconditionerFactory().setAssembledMatrix(A);
        return ParentType::solve(A, x, b);
    }

    /*!
     * \brief Set the matrix A of the linear system Ax = b for reuse
     * \note A is the matrix as assembled; in a parallel run, the solver prepares it (sums its border rows)
     */
    void setMatrix(std::shared_ptr<Matrix> A)
    {
        this->preconditionerFactory().setAssembledMatrix(*A);
        ParentType::setMatrix(std::move(A));
    }

    /*!
     * \brief Set the matrix A of the linear system Ax = b for reuse
     * \note A is the matrix as assembled; the client has to take care of the lifetime management of A
     */
    void setMatrix(Matrix& A)
    { setMatrix(Dune::stackobject_to_shared_ptr(A)); }

    /*!
     * \brief Update the solver after the grid and its dof distribution changed,
     *        e.g. after grid adaption or dynamic load balancing
     * \param gridGeometry the updated grid geometry, which defines the new prolongation
     */
    template<class GridGeometry>
    void updateAfterGridAdaption(const GridGeometry& gridGeometry)
    {
        ParentType::updateAfterGridAdaption(gridGeometry.gridView(), gridGeometry.dofMapper());
        this->setPreconditionerFactory(Factory(gridGeometry));
    }
};

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A conjugate gradient solver with the p-multigrid preconditioner "PQk -> P1 -> AMG"
 *
 * For symmetric positive definite matrices of higher-order Lagrange finite element discretizations
 * (PQ2, PQ3). Preconditioner: PMultigridPreconditioner (SSOR smoothing, Galerkin P1/Q1 coarse level,
 * one algebraic multigrid cycle), see PMultigridPreconditionerFactory for the parameters.
 * \note Constraints imposed by replacing matrix rows (e.g. strongly enforced Dirichlet values) make the matrix
 *       nonsymmetric; eliminating the constrained columns (symmetrizeConstraints()) restores the symmetry.
 * \note In parallel, the smoother takes the corrections of border dofs from the owning process, which makes the
 *       preconditioner slightly nonsymmetric (see PMultigridPreconditionerFactory).
 * \note With the ILU smoother, the preconditioner need not be positive definite (see PMultigridPreconditionerFactory).
 */
template<class LSTraits, class LATraits>
using PMultigridCGIstlSolver =
    Detail::PMultigridIstlSolver<LSTraits, LATraits, Dune::CGSolver<typename LATraits::Vector>>;

/*!
 * \ingroup Linear
 * \brief A flexible conjugate gradient solver (Dune::GeneralizedPCGSolver) with the p-multigrid preconditioner
 *
 * For symmetric positive definite matrices like PMultigridCGIstlSolver. Each search direction is orthogonalized
 * against the previous ones (at most `LinearSolver.Restart`), which keeps the method robust for preconditioners that
 * are not exactly symmetric, such as the p-multigrid preconditioner in parallel, at the cost of the stored directions.
 */
template<class LSTraits, class LATraits>
using PMultigridGeneralizedPCGIstlSolver =
    Detail::PMultigridIstlSolver<LSTraits, LATraits, Dune::GeneralizedPCGSolver<typename LATraits::Vector>>;

/*!
 * \ingroup Linear
 * \brief A BiCGSTAB solver with the p-multigrid preconditioner "PQk -> P1 -> AMG"
 *
 * For the nonsymmetric matrices of hybrid control-volume finite element discretizations (PQ2, PQ3).
 * Preconditioner: PMultigridPreconditioner, see PMultigridPreconditionerFactory for the parameters
 * (for nonsymmetric matrices the coarse AMG may need `LinearSolver.Preconditioner.AmgCriterionSymmetric = false`).
 */
template<class LSTraits, class LATraits>
using PMultigridBiCGSTABIstlSolver =
    Detail::PMultigridIstlSolver<LSTraits, LATraits, Dune::BiCGSTABSolver<typename LATraits::Vector>>;

} // end namespace Dumux

#endif
