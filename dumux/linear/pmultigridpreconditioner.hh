// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Two-level p-multigrid preconditioner for higher-order Lagrange discretizations
 */
#ifndef DUMUX_LINEAR_PMULTIGRID_PRECONDITIONER_HH
#define DUMUX_LINEAR_PMULTIGRID_PRECONDITIONER_HH

#include <bitset>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parametertree.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#if HAVE_MPI
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/schwarz.hh>
#endif

#include <dumux/common/multimapperview.hh>
#include <dumux/linear/lagrangeprolongation.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/pmultigridsmoothers.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Two-level p-multigrid preconditioner "PQk -> P1 -> AMG"
 *
 * The fine level is the operator \f$A\f$ of a higher-order Lagrange discretization. The coarse
 * level is the P1/Q1 space on the same grid (one dof per grid vertex), the prolongation
 * \f$P\f$ interpolates P1/Q1 functions in the fine dofs (see makeLagrangeP1Prolongation()), and the
 * coarse problem with the Galerkin matrix \f$A_c = P^T A P\f$ is solved approximately, usually by one
 * algebraic multigrid cycle. One application \f$v = M^{-1} d\f$ performs
 *  1. smoothing of \f$A v = d\f$ starting from \f$v = 0\f$,
 *  2. restriction of the residual \f$d - A v\f$ with \f$P^T\f$,
 *  3. the coarse solve,
 *  4. prolongation of the coarse correction and addition to \f$v\f$,
 *  5. smoothing of the residual equation, whose correction is added to \f$v\f$.
 * For a symmetric matrix, a symmetric smoother (Jacobi, SSOR, ILU, an exact solve) and a symmetric coarse solver,
 * the preconditioner is symmetric. In parallel, the smoother of PMultigridPreconditionerFactory takes the corrections of border dofs from the
 * owning process (restricted additive Schwarz), which makes the preconditioner slightly nonsymmetric.
 *
 * The parts are passed to the constructor, which makes the class independent of the decomposition:
 * for a parallel operator, the smoother and the coarse solver are parallel preconditioners, one function
 * brings the residual into the representation the restriction needs and another one the restricted defect into
 * that of the coarse solver. PMultigridPreconditionerFactory constructs the parts for sequential, nonoverlapping
 * and overlapping parallel operators.
 *
 * \tparam X the type of the update vector
 * \tparam Y the type of the defect vector
 * \tparam CoarseVector the vector type of the coarse level
 * \tparam ProlongationMatrix the matrix type of the prolongation (fine rows, coarse columns)
 */
template<class X, class Y, class CoarseVector, class ProlongationMatrix>
class PMultigridPreconditioner : public Dune::Preconditioner<X, Y>
{
public:
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;
    using FineOperator = Dune::LinearOperator<X, Y>;
    using CoarseSolver = Dune::Preconditioner<CoarseVector, CoarseVector>;

    /*!
     * \brief Constructor
     * \param fineOperator the operator of the fine level (residual computation)
     * \param smoother the smoother of the fine level, applied from a zero initial guess
     * \param prolongation the prolongation from the coarse to the fine level
     * \param coarseSolver the (approximate) solver of the coarse level
     * \param prepareResidual applied to the residual before the restriction (empty: nothing)
     * \param makeCoarseDefectConsistent applied to the restricted defect (empty: nothing)
     * \param coarseLevelData kept alive as long as the preconditioner (e.g. the coarse matrix and communication)
     */
    PMultigridPreconditioner(std::shared_ptr<const FineOperator> fineOperator,
                             std::shared_ptr<Dune::Preconditioner<X, Y>> smoother,
                             std::shared_ptr<const ProlongationMatrix> prolongation,
                             std::shared_ptr<CoarseSolver> coarseSolver,
                             std::function<void(Y&)> prepareResidual = {},
                             std::function<void(CoarseVector&)> makeCoarseDefectConsistent = {},
                             std::shared_ptr<const void> coarseLevelData = {})
    : coarseLevelData_(std::move(coarseLevelData))
    , makeCoarseDefectConsistent_(std::move(makeCoarseDefectConsistent))
    , fineOperator_(std::move(fineOperator))
    , smoother_(std::move(smoother))
    , prolongation_(std::move(prolongation))
    , coarseSolver_(std::move(coarseSolver))
    , prepareResidual_(std::move(prepareResidual))
    , residual_(prolongation_->N())
    , correction_(prolongation_->N())
    , coarseCorrection_(prolongation_->M())
    , coarseDefect_(prolongation_->M())
    {}

    //! Prepare the smoother and the coarse solver
    void pre(X& x, Y& b) override
    {
        smoother_->pre(x, b);
        prepareCoarseSolver_();
    }

    /*!
     * \brief Apply the preconditioner
     * \param v the update (output)
     * \param d the defect
     */
    void apply(X& v, const Y& d) override
    {
        prepareCoarseSolver_();

        v = 0.0;
        smoother_->apply(v, d);

        computeResidual_(v, d);
        restrictResidual_();
        coarseCorrection_ = 0.0;
        coarseSolver_->apply(coarseCorrection_, coarseDefect_);
        prolongation_->umv(coarseCorrection_, v);

        computeResidual_(v, d);
        correction_ = 0.0;
        smoother_->apply(correction_, residual_);
        v += correction_;
    }

    //! Clean up
    void post(X& x) override
    {
        smoother_->post(x);
        if (coarsePrepared_)
            coarseSolver_->post(coarseCorrection_);
        coarsePrepared_ = false;
    }

    //! The category of the fine operator
    Dune::SolverCategory::Category category() const override
    { return fineOperator_->category(); }

private:
    void prepareCoarseSolver_()
    {
        if (coarsePrepared_)
            return;
        coarseCorrection_ = 0.0;
        coarseDefect_ = 0.0;
        coarseSolver_->pre(coarseCorrection_, coarseDefect_);
        coarsePrepared_ = true;
    }

    void computeResidual_(const X& v, const Y& d)
    {
        residual_ = d;
        fineOperator_->applyscaleadd(-1.0, v, residual_);
    }

    void restrictResidual_()
    {
        if (prepareResidual_)
            prepareResidual_(residual_);

        coarseDefect_ = 0.0;
        prolongation_->mtv(residual_, coarseDefect_);
        if (makeCoarseDefectConsistent_)
            makeCoarseDefectConsistent_(coarseDefect_);
    }

    // declared first, destroyed last: the coarse solver refers to the coarse level data
    std::shared_ptr<const void> coarseLevelData_;
    std::function<void(CoarseVector&)> makeCoarseDefectConsistent_;
    std::shared_ptr<const FineOperator> fineOperator_;
    std::shared_ptr<Dune::Preconditioner<X, Y>> smoother_;
    std::shared_ptr<const ProlongationMatrix> prolongation_;
    std::shared_ptr<CoarseSolver> coarseSolver_;
    std::function<void(Y&)> prepareResidual_;

    Y residual_;
    X correction_;
    CoarseVector coarseCorrection_, coarseDefect_;
    bool coarsePrepared_ = false;
};

namespace Detail::PMultigrid {

/*!
 * \ingroup Linear
 * \brief The Galerkin product \f$P^T A P\f$, leaving out the rows of \f$A\f$ marked in skipRow
 * \note The blocks of the prolongation are multiples of the identity (see makeLagrangeP1Prolongation()).
 *       The pattern contains the diagonal of every row.
 */
template<class CoarseMatrix, class Matrix, class Prolongation>
CoarseMatrix galerkinProduct(const Matrix& A, const Prolongation& P, const std::vector<bool>& skipRow)
{
    using Scalar = typename Prolongation::field_type;
    using Block = typename CoarseMatrix::block_type;
    const std::size_t numCoarse = P.M();

    // the fine rows contributing to each coarse row, with their weights (the transposed prolongation)
    std::vector<std::vector<std::pair<std::size_t, Scalar>>> fineRows(numCoarse);
    for (auto row = P.begin(); row != P.end(); ++row)
        if (skipRow.empty() || !skipRow[row.index()])
            for (auto col = row->begin(); col != row->end(); ++col)
                fineRows[col.index()].emplace_back(row.index(), (*col)[0][0]);

    // pattern (Gustavson's algorithm with one marker per coarse column)
    std::vector<std::vector<std::size_t>> pattern(numCoarse);
    std::vector<std::size_t> marker(numCoarse, numCoarse);
    std::size_t numNonZeros = 0;
    for (std::size_t I = 0; I < numCoarse; ++I)
    {
        marker[I] = I;
        pattern[I].push_back(I);
        for (const auto& [i, weightI] : fineRows[I])
            for (auto a = A[i].begin(); a != A[i].end(); ++a)
                for (auto p = P[a.index()].begin(); p != P[a.index()].end(); ++p)
                    if (marker[p.index()] != I)
                    {
                        marker[p.index()] = I;
                        pattern[I].push_back(p.index());
                    }
        numNonZeros += pattern[I].size();
    }

    CoarseMatrix Ac(numCoarse, numCoarse, numNonZeros, CoarseMatrix::row_wise);
    for (auto row = Ac.createbegin(); row != Ac.createend(); ++row)
        for (const auto J : pattern[row.index()])
            row.insert(J);

    std::vector<Block> accumulator(numCoarse, Block(0.0));
    for (std::size_t I = 0; I < numCoarse; ++I)
    {
        for (const auto& [i, weightI] : fineRows[I])
            for (auto a = A[i].begin(); a != A[i].end(); ++a)
                for (auto p = P[a.index()].begin(); p != P[a.index()].end(); ++p)
                    accumulator[p.index()].axpy(weightI*(*p)[0][0], *a);

        for (const auto J : pattern[I])
        {
            Ac[I][J] = accumulator[J];
            accumulator[J] = 0.0;
        }
    }

    return Ac;
}

//! Linear solver traits of the P1/Q1 space on the grid vertices, for the parallel helpers
template<class GV>
struct VertexLinearSolverTraits
{
    using GridView = GV;
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
    static constexpr int dofCodim = GV::dimension;
    static constexpr bool canCommunicate = Dune::Capabilities::canCommunicate<typename GV::Grid, GV::dimension>::v;
};

/*!
 * \brief Marks the dofs on entities that are neither interior nor border entities
 * \note The same marking as ParallelISTLHelper::isGhost, available without the parallel helper of a solver
 */
template<class GridGeometry>
std::vector<bool> ghostDofs(const GridGeometry& gridGeometry)
{
    constexpr int dim = GridGeometry::GridView::dimension;
    std::vector<bool> isGhost(gridGeometry.numDofs(), false);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        Dune::Hybrid::forEach(std::make_integer_sequence<int, dim+1>{}, [&](auto codim)
        {
            for (int i = 0; i < static_cast<int>(element.subEntities(codim)); ++i)
            {
                const auto entity = element.template subEntity<codim>(i);
                if (entity.partitionType() != Dune::InteriorEntity && entity.partitionType() != Dune::BorderEntity)
                    for (const auto dofIdx : asMultiMapper(gridGeometry.dofMapper()).indices(entity))
                        isGhost[dofIdx] = true;
            }
        });
    }
    return isGhost;
}

} // end namespace Detail::PMultigrid

/*!
 * \ingroup Linear
 * \brief Constructs the p-multigrid preconditioner (PMultigridPreconditioner) for a linear operator
 *
 * The factory holds the prolongation from the P1/Q1 space on the grid vertices to the dofs of a
 * higher-order Lagrange discretization. It is called like the preconditioner factories of dune-istl
 * with the linear operator and a parameter tree, and returns the preconditioner for
 *  - sequential operators: the fine smoother described below and an AMG with SSOR smoother on the Galerkin coarse
 *    matrix,
 *  - nonoverlapping parallel operators (Dune::NonoverlappingSchwarzOperator on a matrix whose border rows are
 *    summed, see prepareMatrixParallel()): the smoother described below and one coarse problem for all processes
 *    with a parallel AMG. The coarse matrix is
 *    computed from the matrix as assembled (setAssembledMatrix()), whose rows add up over the processes:
 *    the products \f$P^T A P\f$ of the processes, without the rows of ghost dofs, add up to the global
 *    Galerkin matrix and are summed at the border vertices like the fine matrix,
 *  - overlapping parallel operators (Dune::OverlappingSchwarzOperator, e.g. YaspGrid with overlap): the same with
 *    the coarse matrix from the local matrix, whose rows are complete for the fine dofs in the support of the basis
 *    functions of the owned vertices.
 *
 * The fine smoother is SSOR, ILU(n) or an exact LU solve. In a parallel run, it acts on the subdomain of each process
 * including dofs the process does not own, with their defect and correction taken from the owning process (restricted
 * additive Schwarz), which makes the preconditioner slightly nonsymmetric. For an overlapping decomposition the
 * subdomain contains the overlap. For a nonoverlapping one it contains the ghost dofs, whose rows the owning processes
 * send, and the dofs of all elements containing a border vertex of the process (VertexPatchOverlap), so that the rows
 * of the owned dofs are complete. In both cases the rows of the dofs not owned lack the couplings to dofs outside the
 * subdomain, whose absolute values enlarge the singular values of their diagonal blocks (enlargeIncompleteRows()),
 * which keeps the subdomain problem well posed, also for a subdomain without Dirichlet constraints; matrix blocks with
 * more than three rows need LAPACK for this. For a symmetric matrix, all three smoothers keep the sequential
 * preconditioner symmetric, but ILU of a matrix that is not an M-matrix (e.g. of PQ2/PQ3 elements) can have small or
 * negative pivots, so that the preconditioner of a symmetric positive definite matrix need not be positive definite,
 * as the conjugate gradient method requires. ILU and LU handle rows with small or negative diagonal entries, e.g. the
 * control-volume rows of the vertex dofs of hybrid PQ3 schemes, where SSOR smoothing needs strong damping or stalls;
 * in parallel, ILU(0) suffices for these schemes on quadrilaterals, while on simplices ILU(1) or LU is needed.
 *
 * Parameters (input file key; key of the `preconditioner` sub-tree of the linear solver parameters):
 *
 * Input file key                                  | Tree key           | Meaning                                    | Default
 * ------------------------------------------------|--------------------|--------------------------------------------|--------
 * LinearSolver.Preconditioner.PMultigridSmoother  | pmultigridSmoother | fine smoother, `ssor`, `ilu` or `lu` (UMFPack) | ssor
 * LinearSolver.Preconditioner.Iterations          | iterations         | SSOR sweeps of the fine smoother before and after the coarse correction | 1
 * LinearSolver.Preconditioner.ILUOrder            | n                  | fill-in level of the ILU smoother          | 0
 * LinearSolver.Preconditioner.Relaxation          | relaxation         | relaxation factor of the SSOR and ILU smoothers | 1.0
 * LinearSolver.Preconditioner.Amg*                | (AMG keys)         | the coarse AMG, e.g. `AmgSmootherIterations`, `AmgMaxLevel`, `AmgCoarsenTarget`, `AmgCriterionSymmetric` (`false` for nonsymmetric matrices) | see dune-istl
 *
 * \tparam GridView the grid view of the discretization
 * \tparam Matrix the (fine) matrix type, a Dune::BCRSMatrix with square Dune::FieldMatrix blocks
 */
template<class GridView, class Matrix>
class PMultigridPreconditionerFactory
{
    using Block = typename Matrix::block_type;
    using Scalar = typename Block::field_type;
    static constexpr int blockSize = Block::rows;
    static_assert(Block::rows == Block::cols, "Square matrix blocks required");

    using VertexTraits = Detail::PMultigrid::VertexLinearSolverTraits<GridView>;
    using VertexMapper = typename VertexTraits::DofMapper;
    using FineVector = Dune::BlockVector<Dune::FieldVector<Scalar, blockSize>>;
    using AbsoluteRowSums = Dune::BlockVector<Dune::FieldVector<Scalar, blockSize>>;
    using Overlap = Detail::PMultigrid::VertexPatchOverlap<GridView, Matrix, FineVector>;
    using Extension = typename Overlap::Extension;

public:
    //! The matrix type of the prolongation and of the coarse level
    using CoarseMatrix = Dune::BCRSMatrix<Block>;
    using ProlongationMatrix = CoarseMatrix;
    //! The vector type of the coarse level
    using CoarseVector = Dune::BlockVector<Dune::FieldVector<Scalar, blockSize>>;

    //! A factory without prolongation, to be replaced before use
    PMultigridPreconditionerFactory() = default;

    /*!
     * \brief Factory for the higher-order Lagrange discretization of a grid geometry
     * \param gridGeometry a grid geometry with the local dof interface (FE or hybrid CVFE PQk)
     */
    template<class GridGeometry>
    explicit PMultigridPreconditionerFactory(const GridGeometry& gridGeometry)
    : prolongation_(std::make_shared<ProlongationMatrix>(makeLagrangeP1Prolongation<blockSize, Scalar>(gridGeometry)))
    {
        if (gridGeometry.gridView().comm().size() > 1)
            parallelData_ = std::make_shared<ParallelData>(gridGeometry);
    }

    //! Factory for sequential operators with a given prolongation
    explicit PMultigridPreconditionerFactory(std::shared_ptr<const ProlongationMatrix> prolongation)
    : prolongation_(std::move(prolongation))
    {}

    /*!
     * \brief Set the matrix as assembled, before it is prepared for a parallel solver (collective)
     * \note Only needed (and only used) for parallel runs on nonoverlapping decompositions. The preconditioners
     *       constructed afterwards use this matrix; copies of the factory made before keep the previous one.
     */
    void setAssembledMatrix(const Matrix& A)
    {
        if (parallelData_ && !parallelData_->isOverlapping)
        {
            const auto& isGhost = parallelData_->isGhostFineDof;
            assembledCoarseMatrix_ = std::make_shared<const CoarseMatrix>(
                Detail::PMultigrid::galerkinProduct<CoarseMatrix>(A, *prolongation_, isGhost)
            );

            Matrix localSmootherMatrix(A);
            AbsoluteRowSums completeGhostRowSums;
            parallelData_->makeSmootherMatrix(localSmootherMatrix, completeGhostRowSums);
            auto smootherSubdomain = std::make_shared<Extension>();
            parallelData_->overlap->extend(localSmootherMatrix, *smootherSubdomain);
            Detail::PMultigrid::enlargeIncompleteRows(smootherSubdomain->matrix, completeGhostRowSums, isGhost);
            smootherSubdomain_ = std::move(smootherSubdomain);
        }
    }

    //! The prolongation from the coarse to the fine level
    const ProlongationMatrix& prolongation() const
    { return *prolongation_; }

    /*!
     * \brief Construct the preconditioner for a linear operator
     * \param op the linear operator, a Dune::MatrixAdapter, Dune::NonoverlappingSchwarzOperator or
     *           Dune::OverlappingSchwarzOperator
     * \param config the parameters (see class documentation)
     * \note The first argument (operator information of the dune-istl factory interface) is not used.
     */
    template<class OpInfo, class Op>
    std::shared_ptr<Dune::Preconditioner<typename Op::domain_type, typename Op::range_type>>
    operator()(OpInfo, const std::shared_ptr<Op>& op, const Dune::ParameterTree& config) const
    {
        using M = typename Op::matrix_type;
        using X = typename Op::domain_type;
        using Y = typename Op::range_type;
        static_assert(std::is_same_v<typename M::block_type, Block>, "Operator and factory have different matrix blocks");
        using Preconditioner = PMultigridPreconditioner<X, Y, CoarseVector, ProlongationMatrix>;

        if (!prolongation_)
            DUNE_THROW(Dune::InvalidStateException, "p-multigrid factory without prolongation");

        const M& A = op->getmat();
        if (A.N() != prolongation_->N())
            DUNE_THROW(Dune::InvalidStateException, "Prolongation rows (" << prolongation_->N()
                       << ") and matrix rows (" << A.N() << ") differ in number");

        if (op->category() == Dune::SolverCategory::sequential)
        {
            using CoarseOperator = Dune::MatrixAdapter<CoarseMatrix, CoarseVector, CoarseVector>;
            using CoarseSmoother = Dune::SeqSSOR<CoarseMatrix, CoarseVector, CoarseVector>;
            using CoarseSolver = Dune::Amg::AMG<CoarseOperator, CoarseVector, CoarseSmoother>;

            auto coarseMatrix = std::make_shared<const CoarseMatrix>(
                Detail::PMultigrid::galerkinProduct<CoarseMatrix>(A, *prolongation_, {})
            );
            auto coarseOperator = std::make_shared<const CoarseOperator>(coarseMatrix);
            auto coarseSolver = std::make_shared<CoarseSolver>(coarseOperator, config);
            return std::make_shared<Preconditioner>(op, makeSequentialSmoother_<M, X, Y>(A, config), prolongation_, coarseSolver);
        }

#if HAVE_MPI
        if constexpr (requires { op->getCommunication(); })
        {
            using Comm = std::decay_t<decltype(op->getCommunication())>;
            using Smoother = Detail::PMultigrid::RestrictedSchwarzSmoother<X, Y, Comm>;
            const auto category = op->category();
            const auto& comm = op->getCommunication();
            if (!parallelData_)
                DUNE_THROW(Dune::InvalidStateException, "The parallel p-multigrid preconditioner needs a factory constructed from the grid geometry");
            const auto& data = *parallelData_;

            if (category == Dune::SolverCategory::nonoverlapping)
            {
                using CoarseOperator = Dune::NonoverlappingSchwarzOperator<CoarseMatrix, CoarseVector, CoarseVector, Comm>;
                using CoarseSmoother = Dune::NonoverlappingBlockPreconditioner<Comm, Dune::SeqSSOR<CoarseMatrix, CoarseVector, CoarseVector>>;
                using CoarseSolver = Dune::Amg::AMG<CoarseOperator, CoarseVector, CoarseSmoother, Comm>;

                if (!assembledCoarseMatrix_ || !smootherSubdomain_)
                    DUNE_THROW(Dune::InvalidStateException, "The nonoverlapping p-multigrid preconditioner needs the matrix as assembled"
                               << " (setAssembledMatrix)");

                struct CoarseLevel
                {
                    std::shared_ptr<CoarseMatrix> matrix;
                    std::shared_ptr<Comm> communication;
                    std::shared_ptr<CoarseOperator> op;
                };
                auto coarse = std::make_shared<CoarseLevel>();
                coarse->matrix = std::make_shared<CoarseMatrix>(*assembledCoarseMatrix_);
                data.makeCoarseMatrixConsistent(*coarse->matrix);
                coarse->communication = std::make_shared<Comm>(data.gridView.comm(), category);
                data.coarseHelper->createParallelIndexSet(*coarse->communication);
                coarse->op = std::make_shared<CoarseOperator>(coarse->matrix, *coarse->communication);
                auto coarseSolver = std::make_shared<CoarseSolver>(coarse->op, config, *coarse->communication);

                // each residual entry is restricted by the process owning its dof
                const auto isNonOwned = Detail::PMultigrid::nonOwnedDofs(comm, A.N());
                std::vector<Scalar> ownerWeights(A.N(), 1.0);
                for (std::size_t i = 0; i < A.N(); ++i)
                    if (isNonOwned[i])
                        ownerWeights[i] = 0.0;
                auto prepareResidual = [weights = std::move(ownerWeights)](Y& r)
                {
                    for (std::size_t i = 0; i < r.size(); ++i)
                        r[i] *= weights[i];
                };

                auto makeConsistent = [communication = coarse->communication,
                                       isGhost = data.isGhostCoarseDof](CoarseVector& d)
                {
                    communication->addOwnerCopyToOwnerCopy(d, d);
                    for (std::size_t I = 0; I < d.size(); ++I)
                        if (isGhost[I])
                            d[I] = 0.0;
                };

                static_assert(std::is_same_v<X, FineVector> && std::is_same_v<Y, FineVector>, "Unexpected vector type");
                using OverlapSmoother = Detail::PMultigrid::VertexPatchSchwarzSmoother<X, Y, Comm, Overlap>;
                auto localSmoother = makeSequentialSmoother_<M, X, Y>(smootherSubdomain_->matrix, config);
                auto smoother = std::make_shared<OverlapSmoother>(
                    localSmoother, smootherSubdomain_, data.overlap, comm
                );
                return std::make_shared<Preconditioner>(op, smoother, prolongation_, coarseSolver,
                                                        std::move(prepareResidual), std::move(makeConsistent), coarse);
            }

            if (category == Dune::SolverCategory::overlapping)
            {
                using CoarseOperator = Dune::OverlappingSchwarzOperator<CoarseMatrix, CoarseVector, CoarseVector, Comm>;
                using CoarseSmoother = Dune::BlockPreconditioner<CoarseVector, CoarseVector, Comm, Dune::SeqSSOR<CoarseMatrix, CoarseVector, CoarseVector>>;
                using CoarseSolver = Dune::Amg::AMG<CoarseOperator, CoarseVector, CoarseSmoother, Comm>;

                // the rows of owned dofs are complete in the local matrix of an overlapping decomposition, and so are the
                // rows of the fine dofs in the support of the basis function of an owned vertex
                struct CoarseLevel
                {
                    std::shared_ptr<CoarseMatrix> matrix;
                    std::shared_ptr<Comm> communication;
                    std::shared_ptr<CoarseOperator> op;
                };
                auto coarse = std::make_shared<CoarseLevel>();
                coarse->matrix = std::make_shared<CoarseMatrix>(Detail::PMultigrid::galerkinProduct<CoarseMatrix>(A, *prolongation_, {}));
                coarse->communication = std::make_shared<Comm>(data.gridView.comm(), category);
                data.coarseHelper->createParallelIndexSet(*coarse->communication);
                coarse->op = std::make_shared<CoarseOperator>(coarse->matrix, *coarse->communication);
                auto coarseSolver = std::make_shared<CoarseSolver>(coarse->op, config, *coarse->communication);

                auto prepareResidual = [&comm](Y& r) { comm.copyOwnerToAll(r, r); };
                auto makeConsistent = [communication = coarse->communication](CoarseVector& d) { communication->project(d); };

                // the rows of the dofs not owned lack the couplings to dofs outside the overlap
                auto smootherMatrix = std::make_shared<M>(A);
                AbsoluteRowSums completeRowSums(A.N());
                for (std::size_t i = 0; i < A.N(); ++i)
                    completeRowSums[i] = Detail::PMultigrid::absoluteRowSum(A[i]);
                comm.copyOwnerToAll(completeRowSums, completeRowSums);
                Detail::PMultigrid::enlargeIncompleteRows(*smootherMatrix, completeRowSums,
                                                          Detail::PMultigrid::nonOwnedDofs(comm, A.N()));

                auto localSmoother = makeSequentialSmoother_<M, X, Y>(*smootherMatrix, config);
                auto smoother = std::make_shared<Smoother>(localSmoother, smootherMatrix, comm, category);
                return std::make_shared<Preconditioner>(op, smoother, prolongation_, coarseSolver,
                                                        std::move(prepareResidual), std::move(makeConsistent), coarse);
            }
        }
#endif

        DUNE_THROW(Dune::NotImplemented, "p-multigrid preconditioner for solver category " << op->category());
    }

    /*!
     * \brief Construct the preconditioner for a linear operator
     * \param op the linear operator, a Dune::MatrixAdapter, Dune::NonoverlappingSchwarzOperator or
     *           Dune::OverlappingSchwarzOperator
     * \param config the parameters (see class documentation)
     */
    template<class Op>
    std::shared_ptr<Dune::Preconditioner<typename Op::domain_type, typename Op::range_type>>
    operator()(const std::shared_ptr<Op>& op, const Dune::ParameterTree& config) const
    { return (*this)(0, op, config); }

private:
    template<class M, class X, class Y>
    static std::shared_ptr<Dune::Preconditioner<X, Y>> makeSequentialSmoother_(const M& A, const Dune::ParameterTree& config)
    {
        const auto relaxation = config.get<Scalar>("relaxation", 1.0);
        const auto smootherType = config.get<std::string>("pmultigridSmoother", "ssor");
        if (smootherType == "ssor")
            return std::make_shared<Dune::SeqSSOR<M, X, Y>>(A, config.get<int>("iterations", 1), relaxation);
        else if (smootherType == "ilu")
            return std::make_shared<Dune::SeqILU<M, X, Y>>(A, config.get<int>("n", 0), relaxation);
        else if (smootherType == "lu")
        {
#if HAVE_SUITESPARSE_UMFPACK
            return std::make_shared<Detail::PMultigrid::DirectSolve<M, X, Y>>(A);
#else
            DUNE_THROW(Dune::InvalidStateException, "The p-multigrid smoother lu needs UMFPack (SuiteSparse)");
#endif
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Unknown p-multigrid smoother " << smootherType << " (ssor, ilu, lu)");
    }

    //! The vertex space of the coarse level in a parallel run
    struct ParallelData
    {
        template<class GridGeometry>
        explicit ParallelData(const GridGeometry& gridGeometry)
        : gridView(gridGeometry.gridView())
        , isOverlapping(gridView.overlapSize(0) > 0)
        , vertexMapper(gridView, Dune::mcmgVertexLayout())
        , isGhostFineDof(Detail::PMultigrid::ghostDofs(gridGeometry))
        , isGhostCoarseDof(vertexMapper.size(), false)
        {
            if constexpr (!VertexTraits::canCommunicate)
                DUNE_THROW(Dune::InvalidStateException, "The parallel p-multigrid preconditioner needs a grid communicating vertices");
            else
            {
                constexpr int dim = GridView::dimension;
                for (const auto& vertex : vertices(gridView))
                    if (vertex.partitionType() != Dune::InteriorEntity && vertex.partitionType() != Dune::BorderEntity)
                        isGhostCoarseDof[vertexMapper.index(vertex)] = true;

                coarseHelper = std::make_shared<ParallelISTLHelper<VertexTraits>>(gridView, vertexMapper);
                coarseMatrixHelper = std::make_shared<CoarseMatrixHelper>(
                    gridView, vertexMapper, std::bitset<dim+1>{ 1UL << dim }
                );

                // the extended mode also sends the rows of the dofs on interior entities next to the process border,
                // which arrive at the ghost dofs of the neighbouring processes
                using DofMapper = std::decay_t<decltype(gridGeometry.dofMapper())>;
                using FineMatrixHelper = MultiCodimParallelMatrixHelper<Matrix, GridView, DofMapper, dim+1>;
                std::bitset<dim+1> dofCodims;
                for (int codim = 0; codim <= dim; ++codim)
                    dofCodims[codim] = !gridGeometry.dofMapper().types(codim).empty();
                auto fineHelper = std::make_shared<FineMatrixHelper>(gridView, gridGeometry.dofMapper(), dofCodims, true);
                if (!isOverlapping)
                    overlap = std::make_shared<Overlap>(gridGeometry, vertexMapper);
                makeSmootherMatrix = [fineHelper, dofMapper = &gridGeometry.dofMapper(), gv = gridView, dofCodims,
                                      isGhost = isGhostFineDof](Matrix& A, AbsoluteRowSums& completeGhostRowSums)
                {
                    for (std::size_t i = 0; i < A.N(); ++i)
                        if (isGhost[i])
                            A[i] = 0.0;
                    fineHelper->extendMatrix(A, [](auto){ return false; });
                    fineHelper->sumEntries(A);

                    // the absolute row sums of the complete rows, sent from the owning processes to the ghost dofs
                    completeGhostRowSums.resize(A.N());
                    completeGhostRowSums = 0.0;
                    for (std::size_t i = 0; i < A.N(); ++i)
                        if (!isGhost[i])
                            completeGhostRowSums[i] = Detail::PMultigrid::absoluteRowSum(A[i]);
                    MultiCodimVectorCommDataHandleEqual<DofMapper, AbsoluteRowSums, dim> ownerToGhost(
                        *dofMapper, completeGhostRowSums, dofCodims
                    );
                    gv.communicate(ownerToGhost, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                };
            }
        }

        // ParallelData is shared between copies of the factory and referenced by the helpers
        ParallelData(const ParallelData&) = delete;
        ParallelData& operator=(const ParallelData&) = delete;

        //! Sum the border rows of the coarse matrix over the processes and set identity rows for ghost vertices
        void makeCoarseMatrixConsistent(CoarseMatrix& Ac) const
        {
            coarseMatrixHelper->extendMatrix(Ac, [this](auto idx){ return isGhostCoarseDof[idx]; });
            coarseMatrixHelper->sumEntries(Ac);
            for (std::size_t I = 0; I < Ac.N(); ++I)
            {
                if (isGhostCoarseDof[I])
                {
                    Ac[I] = 0.0;
                    for (int k = 0; k < blockSize; ++k)
                        Ac[I][I][k][k] = 1.0;
                }
            }
        }

        using CoarseMatrixHelper = MultiCodimParallelMatrixHelper<CoarseMatrix, GridView, VertexMapper, GridView::dimension+1>;

        GridView gridView;
        bool isOverlapping;
        VertexMapper vertexMapper;
        std::vector<bool> isGhostFineDof;
        std::vector<bool> isGhostCoarseDof;
        std::shared_ptr<ParallelISTLHelper<VertexTraits>> coarseHelper;
        std::shared_ptr<CoarseMatrixHelper> coarseMatrixHelper;
        //! turns a copy of the matrix as assembled into the local matrix of the smoother (ghost rows from the owners)
        std::function<void(Matrix&, AbsoluteRowSums&)> makeSmootherMatrix;
        std::shared_ptr<const Overlap> overlap;
    };

    std::shared_ptr<const ProlongationMatrix> prolongation_;
    std::shared_ptr<ParallelData> parallelData_;
    // computed from the matrix as assembled, for nonoverlapping parallel operators
    std::shared_ptr<const CoarseMatrix> assembledCoarseMatrix_;
    std::shared_ptr<const Extension> smootherSubdomain_;
};

} // end namespace Dumux

#endif
