// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Krylov solvers with a block-diagonal AMG preconditioner for the embedded coupling, also on distributed grids
 */
#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_BLOCKDIAG_AMG_SOLVER_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_BLOCKDIAG_AMG_SOLVER_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/version.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/solverfactory.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>
#if HAVE_MPI
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/novlpschwarz.hh>
#include <dune/istl/schwarz.hh>
#endif

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/scalarproducts.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/solvercategory.hh>

namespace Dumux::Detail::EmbeddedBlockDiagAMG {

//! The category of a system whose subdomains have the given categories, as for the multi-type scalar product
template<std::size_t n>
Dune::SolverCategory::Category combinedCategory(const std::array<Dune::SolverCategory::Category, n>& categories)
{
    const auto all = [&](auto category)
    { return std::all_of(categories.begin(), categories.end(), [&](auto c){ return c == category; }); };

    if (all(Dune::SolverCategory::sequential))
        return Dune::SolverCategory::sequential;
    else if (all(Dune::SolverCategory::nonoverlapping))
        return Dune::SolverCategory::nonoverlapping;
    return Dune::SolverCategory::overlapping;
}

#if HAVE_MPI
using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
#else
struct Comm {};
#endif

/*!
 * \brief The communication of the subdomains of a multi-type vector
 *
 * A consistent vector has on every process the value of the owner for each of its degrees of freedom. The result
 * of the operator has the representation of the operator of a single subdomain: on an overlapping decomposition the
 * product of the rows of the process with the copies set to zero, on a non-overlapping one the sum of the partial
 * products of all processes holding the degree of freedom.
 */
template<std::size_t n>
class SubdomainCommunication
{
public:
    SubdomainCommunication(const std::array<std::shared_ptr<const Comm>, n>& comms,
                           const std::array<Dune::SolverCategory::Category, n>& categories)
    : comms_(comms), categories_(categories)
    {}

    //! Copy the values of the owners to all other processes
    template<class Vector>
    void makeConsistent(Vector& v) const
    {
#if HAVE_MPI
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<n>{}), [&](auto i)
        {
            if (categories_[i] != Dune::SolverCategory::sequential)
                comms_[i]->copyOwnerToAll(v[i], v[i]);
        });
#endif
    }

    //! Bring the partial products of the rows of every process into the representation of the subdomain operator
    template<class Vector>
    void completeProduct(Vector& y) const
    {
#if HAVE_MPI
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<n>{}), [&](auto i)
        {
            if (categories_[i] == Dune::SolverCategory::overlapping)
                comms_[i]->project(y[i]);
            else if (categories_[i] == Dune::SolverCategory::nonoverlapping)
            {
                for (const auto& pair : comms_[i]->indexSet())
                    if (pair.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::overlap)
                        y[i][pair.local().local()] = 0.0;
                comms_[i]->addOwnerCopyToOwnerCopy(y[i], y[i]);
            }
        });
#endif
    }

private:
    std::array<std::shared_ptr<const Comm>, n> comms_;
    std::array<Dune::SolverCategory::Category, n> categories_;
};

/*!
 * \brief The linear operator of the coupled system: the local matrix and the Jacobian entries of the coupling whose
 *        column is a degree of freedom of another process
 */
template<class Matrix, class Vector, class CouplingManager, std::size_t n>
class CoupledOperator : public Dune::LinearOperator<Vector, Vector>
{
public:
    using matrix_type = Matrix;
    using domain_type = Vector;
    using range_type = Vector;
    using field_type = typename Vector::field_type;

    CoupledOperator(const Matrix& A, CouplingManager& couplingManager,
                    const SubdomainCommunication<n>& communication,
                    Dune::SolverCategory::Category category)
    : A_(A), couplingManager_(couplingManager), communication_(communication), category_(category)
    {}

    void apply(const Vector& x, Vector& y) const override
    {
        Vector xc(x);
        communication_.makeConsistent(xc);
        y = 0.0;
        A_.umv(xc, y);
        couplingManager_.addRemoteJacobianProduct(xc, y);
        communication_.completeProduct(y);
    }

    void applyscaleadd(field_type alpha, const Vector& x, Vector& y) const override
    {
        Vector z(y);
        apply(x, z);
        y.axpy(alpha, z);
    }

    Dune::SolverCategory::Category category() const override
    { return category_; }

private:
    const Matrix& A_;
    CouplingManager& couplingManager_;
    const SubdomainCommunication<n>& communication_;
    Dune::SolverCategory::Category category_;
};

//! A preconditioner per subdomain applied to the subdomain blocks of the vectors
template<class Vector, class Preconditioners, std::size_t n>
class BlockDiagonalPreconditioner : public Dune::Preconditioner<Vector, Vector>
{
public:
    BlockDiagonalPreconditioner(const Preconditioners& preconditioners, Dune::SolverCategory::Category category)
    : preconditioners_(preconditioners), category_(category)
    {}

    void pre(Vector& x, Vector& b) override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<n>{}), [&](auto i){ std::get<i>(preconditioners_)->pre(x[i], b[i]); });
    }

    void apply(Vector& v, const Vector& d) override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<n>{}), [&](auto i){ std::get<i>(preconditioners_)->apply(v[i], d[i]); });
    }

    void post(Vector& x) override
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<n>{}), [&](auto i){ std::get<i>(preconditioners_)->post(x[i]); });
    }

    Dune::SolverCategory::Category category() const override
    { return category_; }

private:
    Preconditioners preconditioners_;
    Dune::SolverCategory::Category category_;
};

//! An AMG preconditioner for a linear operator of a single subdomain, as in the AMG solvers of a single domain
template<class LinearOperator>
auto makeAmg(std::shared_ptr<LinearOperator> op, const Dune::ParameterTree& params)
{
    using Prec = Dune::Preconditioner<typename LinearOperator::domain_type, typename LinearOperator::range_type>;
#if DUNE_VERSION_GTE(DUNE_ISTL,2,11)
    using OpTraits = Dune::OperatorTraits<LinearOperator>;
    std::shared_ptr<Prec> prec = Dune::AMGCreator{}(OpTraits{}, op, params);
#if HAVE_MPI
    if constexpr (OpTraits::isParallel)
    {
        using C = typename OpTraits::comm_type;
        const C& comm = OpTraits::getCommOrThrow(op);
        if (op->category() == Dune::SolverCategory::overlapping && prec->category() == Dune::SolverCategory::sequential)
            prec = std::make_shared<Dune::BlockPreconditioner<typename OpTraits::domain_type, typename OpTraits::range_type, C>>(prec, comm);
        else if (op->category() == Dune::SolverCategory::nonoverlapping && prec->category() == Dune::SolverCategory::sequential)
            prec = std::make_shared<Dune::NonoverlappingBlockPreconditioner<C, Prec>>(prec, comm);
    }
#endif
#else
    using TL = Dune::TypeList<typename LinearOperator::matrix_type, typename LinearOperator::domain_type, typename LinearOperator::range_type>;
    std::shared_ptr<Prec> prec = Dune::AMGCreator{}(TL{}, op, params);
#if HAVE_MPI
    if (prec->category() != op->category() && prec->category() == Dune::SolverCategory::sequential)
        prec = Dune::wrapPreconditioner4Parallel(prec, op);
#endif
#endif
    return prec;
}

} // end namespace Dumux::Detail::EmbeddedBlockDiagAMG

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief A Krylov solver with a block-diagonal AMG preconditioner for the embedded coupling of two domains
 *
 * The operator is the assembled Jacobian plus the Jacobian entries of the coupling whose column is a degree of
 * freedom of another process (`addRemoteJacobianProduct` of the coupling manager). The preconditioner applies an AMG
 * to each diagonal block with the communication of its domain, set up as in the AMG solvers of a single domain, so
 * it neglects the off-diagonal blocks and is meant for weak coupling. On one process the same solver runs
 * sequentially. The domains may have different decompositions, e.g. a non-overlapping one for a box discretization
 * and an overlapping one for a cell-centered discretization with overlap.
 *
 * Parameters (group `LinearSolver`) as for the other dune-istl solvers; the Krylov solver reads its parameters from
 * the tree of the first domain, the AMG of each domain from the tree of its domain.
 *
 * \tparam CouplingManager the embedded coupling manager
 * \tparam InverseOperator the dune-istl Krylov solver, e.g. `Dune::BiCGSTABSolver` or `Dune::RestartedGMResSolver`
 */
template<class CouplingManager, template<class...> class InverseOperator>
class EmbeddedBlockDiagAMGSolver : public LinearSolver
{
    using MDTraits = typename CouplingManager::MultiDomainTraits;
    static constexpr std::size_t numDomains = MDTraits::numSubDomains;
    static_assert(numDomains == 2, "The embedded coupling has two domains");

    template<std::size_t i>
    using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i>
    using Traits = Dumux::LinearSolverTraits<GridGeometry<i>>;

    using Comm = Detail::EmbeddedBlockDiagAMG::Comm;
    using Communication = Detail::EmbeddedBlockDiagAMG::SubdomainCommunication<numDomains>;
    using GridGeometries = std::tuple<std::shared_ptr<const GridGeometry<0>>, std::shared_ptr<const GridGeometry<1>>>;

public:
    using Matrix = typename MDTraits::JacobianMatrix;
    using Vector = typename MDTraits::ResidualVector;
    using Scalar = typename Vector::field_type;

private:
    template<std::size_t i>
    using BlockMatrix = std::decay_t<decltype(std::declval<Matrix>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;
    template<std::size_t i>
    using BlockVector = std::decay_t<decltype(std::declval<Vector>()[Dune::index_constant<i>{}])>;
    template<std::size_t i>
    using BlockPreconditioner = Dune::Preconditioner<BlockVector<i>, BlockVector<i>>;

    template<std::size_t i>
    struct Domain
    {
        Dune::SolverCategory::Category category = Dune::SolverCategory::sequential;
        Dune::ParameterTree params;
#if HAVE_MPI
        std::shared_ptr<Comm> comm;
        std::shared_ptr<ParallelISTLHelper<Traits<i>>> helper;
#endif
        std::shared_ptr<BlockMatrix<i>> matrix;
        std::shared_ptr<Dune::AssembledLinearOperator<BlockMatrix<i>, BlockVector<i>, BlockVector<i>>> op;
        std::shared_ptr<BlockPreconditioner<i>> preconditioner;
    };

public:
    /*!
     * \brief Construct the solver
     * \param gridGeometries the grid geometries of the two domains
     * \param couplingManager the coupling manager, which provides the entries whose column is on another process
     * \param paramGroup the parameter group
     */
    template<class GG0, class GG1>
    EmbeddedBlockDiagAMGSolver(const std::tuple<std::shared_ptr<GG0>, std::shared_ptr<GG1>>& gridGeometries,
                               std::shared_ptr<CouplingManager> couplingManager,
                               const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , gridGeometries_(std::get<0>(gridGeometries), std::get<1>(gridGeometries))
    , couplingManager_(std::move(couplingManager))
    {
        std::array<std::shared_ptr<const Comm>, numDomains> comms;
        std::array<Dune::SolverCategory::Category, numDomains> categories;
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<numDomains>{}), [&](auto i)
        {
            auto& domain = std::get<i>(domains_);
            const auto& gridGeometry = *std::get<i>(gridGeometries_);
            domain.params = LinearSolverParameters<Traits<i>>::createParameterTree(paramGroup);
#if HAVE_MPI
            domain.category = Detail::solverCategory<Traits<i>>(gridGeometry.gridView());
            if constexpr (Traits<i>::canCommunicate)
            {
                if (domain.category != Dune::SolverCategory::sequential)
                {
                    domain.comm = std::make_shared<Comm>(gridGeometry.gridView().comm(), domain.category);
                    domain.helper = std::make_shared<ParallelISTLHelper<Traits<i>>>(gridGeometry.gridView(), gridGeometry.dofMapper());
                    domain.helper->createParallelIndexSet(*domain.comm);
                }
            }
            comms[i] = domain.comm;
            if (gridGeometry.gridView().comm().rank() != 0)
                LinearSolverParameters<Traits<i>>::disableVerbosity(domain.params);
#endif
            categories[i] = domain.category;
        });

        category_ = Detail::EmbeddedBlockDiagAMG::combinedCategory(categories);
        communication_ = std::make_unique<Communication>(comms, categories);
        params_ = std::get<0>(domains_).params;

#if HAVE_MPI
        if (category_ != Dune::SolverCategory::sequential)
            scalarProduct_ = std::make_shared<ParallelMultiTypeScalarProduct<Vector, Comm>>(comms);
        else
#endif
            scalarProduct_ = std::make_shared<Dune::SeqScalarProduct<Vector>>();
    }

    /*!
     * \brief Solve the linear system Ax = b
     * \note Collective on distributed grids. The diagonal blocks of A are used as assembled; b is not changed.
     */
    Detail::IstlSolverResult solve(Matrix& A, Vector& x, const Vector& b)
    {
        Vector rhs(b);
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::index_constant<numDomains>{}), [&](auto i)
        {
            setupPreconditioner_(A, rhs, i);
        });

        using Operator = Detail::EmbeddedBlockDiagAMG::CoupledOperator<Matrix, Vector, CouplingManager, numDomains>;
        auto op = std::make_shared<Operator>(A, *couplingManager_, *communication_, category_);

        using Preconditioners = std::tuple<std::shared_ptr<BlockPreconditioner<0>>, std::shared_ptr<BlockPreconditioner<1>>>;
        using Preconditioner = Detail::EmbeddedBlockDiagAMG::BlockDiagonalPreconditioner<Vector, Preconditioners, numDomains>;
        auto preconditioner = std::make_shared<Preconditioner>(
            Preconditioners{std::get<0>(domains_).preconditioner, std::get<1>(domains_).preconditioner}, category_
        );

        InverseOperator<Vector> solver(op, scalarProduct_, preconditioner, params_);
        Dune::InverseOperatorResult result;
        solver.apply(x, rhs, result);
        communication_->makeConsistent(x);
        return result;
    }

    /*!
     * \brief The 2-norm of a vector with the representation of the residual
     * \note Collective on distributed grids. Every degree of freedom counts once; on a non-overlapping decomposition
     *       the partial values of the processes holding it are summed first.
     */
    Scalar norm(const Vector& x) const
    {
#if HAVE_MPI
        if (category_ != Dune::SolverCategory::sequential)
        {
            Vector y(x);
            using namespace Dune::Hybrid;
            forEach(integralRange(Dune::index_constant<numDomains>{}), [&](auto i)
            {
                const auto& domain = std::get<i>(domains_);
                if (domain.category == Dune::SolverCategory::nonoverlapping)
                {
                    using PT = typename Traits<i>::template ParallelNonoverlapping<BlockMatrix<i>, BlockVector<i>>;
                    prepareVectorParallel<Traits<i>, PT>(y[i], *domain.helper);
                }
            });
            return scalarProduct_->norm(y);
        }
#endif
        return scalarProduct_->norm(x);
    }

    std::string name() const
    { return "block-diagonal AMG preconditioned Krylov solver for the embedded coupling"; }

private:
    //! The AMG of domain i, built from the prepared diagonal block, and the right-hand side of domain i
    template<std::size_t i>
    void setupPreconditioner_(Matrix& A, Vector& rhs, Dune::index_constant<i> id)
    {
        auto& domain = std::get<i>(domains_);
        auto& block = A[id][id];
        const auto& params = domain.params.sub("preconditioner");

        if (domain.category == Dune::SolverCategory::sequential)
        {
            domain.matrix = Dune::stackobject_to_shared_ptr(block);
            using Op = typename Traits<i>::template Sequential<BlockMatrix<i>, BlockVector<i>>::LinearOperator;
            auto op = std::make_shared<Op>(domain.matrix);
            domain.preconditioner = Detail::EmbeddedBlockDiagAMG::makeAmg(op, params);
            domain.op = op;
        }
#if HAVE_MPI
        else if constexpr (Traits<i>::canCommunicate)
        {
            if (domain.category == Dune::SolverCategory::nonoverlapping)
            {
                // the operator applies the block as assembled, the AMG needs the summed border rows
                domain.matrix = std::make_shared<BlockMatrix<i>>(block);
                using PT = typename Traits<i>::template ParallelNonoverlapping<BlockMatrix<i>, BlockVector<i>>;
                prepareMatrixParallel<Traits<i>, PT>(*domain.matrix, *domain.helper);
                prepareVectorParallel<Traits<i>, PT>(rhs[id], *domain.helper);
                auto op = std::make_shared<typename PT::LinearOperator>(domain.matrix, *domain.comm);
                domain.preconditioner = Detail::EmbeddedBlockDiagAMG::makeAmg(op, params);
                domain.op = op;
            }
            else
            {
                // the defect has no values at the copies, as the product of the operator
                domain.comm->project(rhs[id]);

                // the rows of the degrees of freedom of other processes are incomplete and can be singular,
                // e.g. for advection out of the local grid, so the AMG gets Dirichlet rows there
                domain.matrix = std::make_shared<BlockMatrix<i>>(block);
                for (const auto& pair : domain.comm->indexSet())
                {
                    if (pair.local().attribute() == Dune::OwnerOverlapCopyAttributeSet::owner)
                        continue;
                    const auto dof = pair.local().local();
                    auto& row = (*domain.matrix)[dof];
                    row = 0.0;
                    for (std::size_t k = 0; k < row[dof].N(); ++k)
                        row[dof][k][k] = 1.0;
                }

                using PT = typename Traits<i>::template ParallelOverlapping<BlockMatrix<i>, BlockVector<i>>;
                auto op = std::make_shared<typename PT::LinearOperator>(domain.matrix, *domain.comm);
                domain.preconditioner = Detail::EmbeddedBlockDiagAMG::makeAmg(op, params);
                domain.op = op;
            }
        }
#endif
    }

    GridGeometries gridGeometries_;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::tuple<Domain<0>, Domain<1>> domains_;
    Dune::SolverCategory::Category category_;
    std::unique_ptr<Communication> communication_;
    Dune::ParameterTree params_;
    std::shared_ptr<Dune::ScalarProduct<Vector>> scalarProduct_;
};

/*!
 * \ingroup EmbeddedCoupling
 * \brief BiCGSTAB with a block-diagonal AMG preconditioner for the embedded coupling
 */
template<class CouplingManager>
using EmbeddedBlockDiagAMGBiCGSTABSolver = EmbeddedBlockDiagAMGSolver<CouplingManager, Dune::BiCGSTABSolver>;

/*!
 * \ingroup EmbeddedCoupling
 * \brief Restarted GMRes with a block-diagonal AMG preconditioner for the embedded coupling
 */
template<class CouplingManager>
using EmbeddedBlockDiagAMGRestartedGMResSolver = EmbeddedBlockDiagAMGSolver<CouplingManager, Dune::RestartedGMResSolver>;

} // end namespace Dumux

#endif
