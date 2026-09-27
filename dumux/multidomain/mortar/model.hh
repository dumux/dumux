// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Class to hold a decomposition and subdomain solvers to compose a mortar coupling model.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_MODEL_HH
#define DUMUX_MULTIDOMAIN_MORTAR_MODEL_HH

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <memory>
#include <numeric>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/concepts/mortartrace_.hh>

#include "couplingmode.hh"
#include "trace.hh"
#include "decomposition.hh"
#include "solverinterface.hh"
#include "projectorinterface.hh"
#include "projectors.hh"

namespace Dumux::Mortar {

// forward declaration
template<typename MortarSolutionVector,
         typename MortarGridGeometry,
         typename... SubDomainGridGeometries>
class ModelFactory;

/*!
 * \ingroup MortarCoupling
 * \brief Holds a decomposition and associated subdomain solvers to compose a mortar coupling
 *        model, following the mortar methods of \cite Boon2022 \cite Boon2023.
 * \note Construct an instance of this class using a `ModelFactory`.
 */
template<typename MortarSolutionVector,
         typename MortarGridGeometry,
         typename... SubDomainGridGeometries>
class Model
{
    using MortarGrid = typename MortarGridGeometry::GridView::Grid;

 public:
    using SolutionVector = MortarSolutionVector;
    using Decomposition = Mortar::Decomposition<MortarGridGeometry, SubDomainGridGeometries...>;
    using SolverVariant = std::variant<std::shared_ptr<SubDomainSolver<MortarSolutionVector, MortarGrid, SubDomainGridGeometries>>...>;

    /*!
     * \brief Impose the given mortar datum on all subdomains: project each mortar's part of
     *        it onto the traces of the subdomains it couples to, in natural mode with the
     *        orientation sign of each side.
     */
    void setMortar(const MortarSolutionVector& x)
    {
        decomposition_.visitMortars([&] (const auto& mortar) {
            const auto mortarId = decomposition_.id(*mortar);
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] (const auto& subDomain) {
                visitSolverFor_(*subDomain, [&] (auto& solver) {
                    auto trace = getProjector_(*mortar, *subDomain).toTrace(
                        extractEntriesFor(*mortar, x)
                    );
                    // a natural (flux) mortar datum is single-valued with respect to the
                    // mortar's reference normal, so the two sides impose it with opposite signs
                    if (couplingMode_ == CouplingMode::natural)
                        trace *= orientation(decomposition_.id(*subDomain), mortarId);
                    solver.setTraceVariables(mortarId, std::move(trace));
                });
            });
        });
    }

    //! Set the mode in which mortar data enters the subdomains, for all subdomain solvers
    void setCouplingMode(CouplingMode mode)
    {
        couplingMode_ = mode;
        for (auto& variant : solvers_)
            std::visit([&] (auto& s) { s->setCouplingMode(mode); }, variant);
    }

    //! Return the mode in which mortar data currently enters the subdomains
    CouplingMode couplingMode() const
    { return couplingMode_; }

    /*!
     * \brief Set whether the subdomains are solved without the external data of their own
     *        problems, which turns a solve of the model into the action of the linear
     *        interface operator on the mortar data.
     */
    void setHomogeneous(bool homogeneous)
    {
        homogeneous_ = homogeneous;
        for (auto& variant : solvers_)
            std::visit([&] (auto& s) { s->setHomogeneous(homogeneous); }, variant);
    }

    //! Return true while the subdomains are solved without the external data of their problems
    bool isHomogeneous() const
    { return homogeneous_; }

    //! Solve all subdomains, each independently, with the mortar datum last imposed
    void solveSubDomains()
    {
        for (auto& variant : solvers_)
            std::visit([] (auto& s) { s->solve(); }, variant);
    }

    /*!
     * \brief Add the conjugate traces of all subdomains, tested against the mortar basis, to
     *        the given residual, which is resized to the number of mortar degrees of freedom.
     */
    void assembleMortarResidual(MortarSolutionVector& residual) const
    {
        residual.resize(numMortarDofs_);
        decomposition_.visitMortars([&] (const auto& mortar) {
            const auto mortarId = decomposition_.id(*mortar);
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] (const auto& subDomain) {
                visitSolverFor_(*subDomain, [&] (const auto& solver) {
                    const auto vars = solver.assembleTraceVariables(mortarId);
                    auto projected = getProjector_(*mortar, *subDomain).fromTrace(vars);
                    if (projected.size() != mortar->numDofs())
                        DUNE_THROW(Dune::InvalidStateException, "Trace does not have the expected number of entries.");
                    // in natural mode the residual is a value jump, the difference rather than
                    // the sum of the two value traces, so the read carries the same signs as
                    // the imposition; the signs cancel in the operator's linear part but not
                    // in the right-hand side
                    if (couplingMode_ == CouplingMode::natural)
                        projected *= orientation(decomposition_.id(*subDomain), mortarId);
                    std::ranges::for_each(projected, [&, i=std::size_t{0}] (const auto& entry) mutable {
                        residual[mortarDofOffsets_[mortarId] + i++] += entry;
                    });
                });
            });
        });
    }

    //! Return the underlying decomposition
    const Decomposition& decomposition() const
    { return decomposition_; }

    //! Return the total number of degrees of freedom on the entire mortar domain
    std::size_t numMortarDofs() const
    { return numMortarDofs_; }

    //! Return the offset of the given mortar's degrees of freedom within the global mortar vector
    std::size_t mortarDofOffset(const MortarGridGeometry& mortar) const
    { return mortarDofOffsets_[decomposition_.id(mortar)]; }

    /*!
     * \brief Orientation sign of a subdomain relative to the reference normal of a mortar.
     *
     * The reference normal points out of the subdomain with sign +1 and into the one with
     * sign -1; the assignment follows the visitation order of the decomposition. Returns 0
     * for an uncoupled pair.
     */
    int orientation(std::size_t subDomainId, std::size_t mortarId) const
    {
        const auto it = orientations_.find(subDomainId);
        if (it == orientations_.end())
            return 0;
        const auto jt = it->second.find(mortarId);
        return jt == it->second.end() ? 0 : jt->second;
    }

    //! Visit each mortar-subdomain coupling with the mortar, the subdomain solver and the subdomain id
    template<typename Visitor>
    void visitCouplings(Visitor&& v) const
    {
        decomposition_.visitMortars([&] (const auto& mortar) {
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] (const auto& sd) {
                visitSolverFor_(*sd, [&] (const auto& solver) {
                    v(*mortar, solver, decomposition_.id(*sd));
                });
            });
        });
    }

    //! Extract the entries of an individual mortar subdomain from the given solution vector
    MortarSolutionVector extractEntriesFor(const MortarGridGeometry& mortar, const MortarSolutionVector& x) const
    {
        if (x.size() != numMortarDofs_)
            DUNE_THROW(Dune::InvalidStateException, "Given vector does not have the expected length");
        const auto mortarId = decomposition_.id(mortar);
        MortarSolutionVector restricted(mortar.numDofs());
        for (std::size_t i = 0; i < mortar.numDofs(); ++i)
            restricted[i] = x[mortarDofOffsets_[mortarId] + i];
        return restricted;
    }

 private:
    friend ModelFactory<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>;

    template<typename ProjectorFactory>
    Model(Decomposition&& decomposition, std::vector<SolverVariant> solvers, ProjectorFactory&& projectorFactory)
    : numMortarDofs_{0}
    , decomposition_{std::move(decomposition)}
    , solvers_{std::move(solvers)}
    {
        if (decomposition_.numberOfSubDomains() != solvers_.size())
            DUNE_THROW(Dune::InvalidStateException, "Number of solvers and subdomains do not match.");

        mortarDofOffsets_.resize(decomposition_.numberOfMortars());
        decomposition_.visitMortars([&] (const auto& mortar) {
            mortarDofOffsets_[decomposition_.id(*mortar)] = mortar->numDofs();
            numMortarDofs_ += mortar->numDofs();

            const auto id = decomposition_.id(*mortar);
            int sideCount = 0;
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] (const auto& sd) {
                if (sideCount > 1)
                    DUNE_THROW(Dune::InvalidStateException, "A mortar with more than two subdomains has no orientation");
                orientations_[decomposition_.id(*sd)][id] = sideCount++ == 0 ? 1 : -1;
            });
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] <typename SD> (const std::shared_ptr<const SD>& sd) {
                decomposition_.visitSubDomainTraceWith(*mortar, *sd,
                    [&] <typename T> (const std::shared_ptr<T>& tracePtr) {
                        if constexpr (std::is_same_v<std::remove_const_t<T>, Trace<SD, MortarGrid>>) {
                            visitSolverFor_(*sd, [&] (auto& solver) {
                                solver.registerMortarTrace(tracePtr, id);
                                std::derived_from<Projector<MortarSolutionVector>> auto p = projectorFactory(
                                    *mortar, solver, *tracePtr
                                );
                                projectors_.emplace_back(std::make_unique<std::remove_cvref_t<decltype(p)>>(std::move(p)));
                                projectorMap_[decomposition_.id(*sd)][decomposition_.id(*mortar)] = projectors_.size() - 1;
                            });
                        } else {
                            DUNE_THROW(Dune::InvalidStateException, "The visited trace is not the trace kind of the subdomain");
                        }
                    }
                );
            });
        });
        std::exclusive_scan(
            mortarDofOffsets_.begin(), mortarDofOffsets_.end(),
            mortarDofOffsets_.begin(), std::size_t{0}
        );
    }

    template<typename GridGeometry, typename Visitor>
        requires(std::disjunction_v<std::is_same<GridGeometry, SubDomainGridGeometries>...>)
    void visitSolverFor_(const GridGeometry& subDomain, Visitor&& v) const
    {
        for (const auto& variant : solvers_)
            if(
                std::visit([&] <typename Solver> (const std::shared_ptr<Solver>& solver) {
                    if constexpr (std::is_same_v<typename Solver::GridGeometry, GridGeometry>)
                        if (solver->gridGeometry().get() == &subDomain)
                        { v(*solver); return true; }
                    return false;
                }, variant)
            )
                return;
        DUNE_THROW(Dune::InvalidStateException, "Could not find solver matching the given subdomain");
    }

    template<typename GridGeometry>
        requires(std::disjunction_v<std::is_same<GridGeometry, SubDomainGridGeometries>...>)
    const auto& getProjector_(const MortarGridGeometry& mortar, const GridGeometry& subDomain) const
    { return *projectors_.at(projectorMap_.at(decomposition_.id(subDomain)).at(decomposition_.id(mortar))); }

    std::size_t numMortarDofs_;
    Decomposition decomposition_;
    std::vector<SolverVariant> solvers_;
    std::vector<std::size_t> mortarDofOffsets_;
    std::vector<std::unique_ptr<Projector<MortarSolutionVector>>> projectors_;
    std::unordered_map<std::size_t, std::unordered_map<std::size_t, std::size_t>> projectorMap_;
    std::unordered_map<std::size_t, std::unordered_map<std::size_t, int>> orientations_;
    CouplingMode couplingMode_ = CouplingMode::essential;
    bool homogeneous_ = false;
};

/*!
 * \ingroup MortarCoupling
 * \brief Factory for constructing a mortar model.
 */
template<typename MortarSolutionVector,
         typename MortarGridGeometry,
         typename... SubDomainGridGeometries>
class ModelFactory
{
    using MortarGrid = typename MortarGridGeometry::GridView::Grid;
    using SolverVariant = typename Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>::SolverVariant;

    template<typename T>
    static constexpr bool isSupportedSolver
        = std::disjunction_v<std::is_same<typename T::GridGeometry, SubDomainGridGeometries>...>
        and std::derived_from<T, SubDomainSolver<MortarSolutionVector, MortarGrid, typename T::GridGeometry>>;

 public:
     //! Insert a mortar domain
    void insertMortar(std::shared_ptr<const MortarGridGeometry> gg)
    { decompositionFactory_.insertMortar(gg); }

    //! Insert a mortar domain and return this factory
    ModelFactory& withMortar(std::shared_ptr<const MortarGridGeometry> gg)
    { insertMortar(gg); return *this; }

    //! Insert a subdomain solver
    template<typename Solver> requires(!std::is_const_v<Solver> and isSupportedSolver<Solver>)
    void insertSubDomain(std::shared_ptr<Solver> solver)
    {
        decompositionFactory_.insertSubDomain(solver->gridGeometry());
        solvers_.emplace_back(std::move(solver));
    }

    //! Insert a subdomain solver and return this factory
    template<typename Solver> requires(!std::is_const_v<Solver> and isSupportedSolver<Solver>)
    ModelFactory& withSubDomain(std::shared_ptr<Solver> gg)
    { insertSubDomain(gg); return *this; }

    /*!
     * \brief Create a model from all inserted mortars & subdomains using default projectors
     *        (uses FVDefaultProjector instances)
     *
     * Each subdomain solver declares the order of the trace data it accepts, and the
     * projector imposing the mortar on its trace is built for that order. The residual
     * trace is piecewise constant in every case.
     */
    Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...> make() const
    {
        return make([] (const auto& mortarGridGeometry, const auto& solver, const auto& trace) {
            return FVDefaultProjector<MortarSolutionVector>{mortarGridGeometry, trace, solver.traceDataOrder()};
        });
    }

    /*!
     * \brief Create a model from all inserted mortars & subdomains with a custom projector
     *        factory, invoked per coupling with the mortar grid geometry, the subdomain
     *        solver and the trace of the subdomain on the mortar, and returning the
     *        Projector for the pair. The data it produces for the trace must have the
     *        layout the solver accepts.
     */
    template<typename F>
    Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...> make(F&& projectorFactory) const
    { return {decompositionFactory_.make(), solvers_, std::forward<F>(projectorFactory)}; }

 private:
    DecompositionFactory<MortarGridGeometry, SubDomainGridGeometries...> decompositionFactory_;
    std::vector<SolverVariant> solvers_;
};

} // end namespace Dumux::Mortar

#endif
