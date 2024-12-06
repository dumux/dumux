// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Class to hold a decomposition and subdomain solvers to compose a mortar coupling model.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_MODEL_HH
#define DUMUX_MULTIDOMAIN_MORTAR_MODEL_HH

#include <numeric>
#include <memory>
#include <utility>
#include <variant>
#include <concepts>

#include <dumux/parallel/parallel_for.hh>

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
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Holds a decomposition and associated subdomain solvers to compose a mortar coupling model.
 * \note Construct an instance of this class using a `ModelFactory`.
 */
template<typename MortarSolutionVector,
         typename MortarGridGeometry,
         typename... SubDomainGridGeometries>
class Model
{
    using MortarGrid = typename MortarGridGeometry::GridView::Grid;

 public:
    using Decomposition = Mortar::Decomposition<MortarGridGeometry, SubDomainGridGeometries...>;
    using SolverVariant = std::variant<std::shared_ptr<SubDomainSolver<MortarSolutionVector, MortarGrid, SubDomainGridGeometries>>...>;

    void setMortar(const MortarSolutionVector& x)
    {
        decomposition_.visitMortars([&] (const auto& mortar) {
            const auto mortarId = decomposition_.id(*mortar);
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] (const auto& subDomain) {
                visitSolverFor_(*subDomain, [&] (auto& solver) {
                    solver.setTraceVariables(mortarId, getProjector_(*mortar, *subDomain).toTrace(
                        extractEntriesFor(*mortar, x)
                    ));
                });
            });
        });
    }

    void solveSubDomains()
    {
        parallelFor(solvers_.size(), [&] (std::size_t i) {
            std::visit([] (auto& s) { s->solve(); }, solvers_.at(i));
        });
    }

    void assembleMortarResidual(MortarSolutionVector& residual) const
    {
        residual.resize(numMortarDofs_);
        decomposition_.visitMortars([&] (const auto& mortar) {
            const auto mortarId = decomposition_.id(*mortar);
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] (const auto& subDomain) {
                visitSolverFor_(*subDomain, [&] (const auto& solver) {
                    const auto vars = solver.assembleTraceVariables(mortarId);
                    const auto projected = getProjector_(*mortar, *subDomain).fromTrace(vars);
                    if (projected.size() != mortar->numDofs())
                        DUNE_THROW(Dune::InvalidStateException, "Trace does not have the expected number of entries.");
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
            decomposition_.visitCoupledSubDomainsOf(*mortar, [&] <typename SD> (const std::shared_ptr<const SD>& sd) {
                decomposition_.visitSubDomainTraceWith(*mortar, *sd,
                    [&] <typename TraceGrid> (const std::shared_ptr<TraceGrid>& tracePtr) {
                        if constexpr (std::is_same_v<SD, typename TraceGrid::DomainGridGeometry>) {
                            visitSolverFor_(*sd, [&] (auto& solver) {
                                solver.registerMortarTrace(tracePtr, id);
                            });
                            std::derived_from<Projector<MortarSolutionVector>> auto p = projectorFactory(
                                *mortar, *sd, *tracePtr, decomposition
                            );
                            projectors_.emplace_back(std::make_unique<std::remove_cvref_t<decltype(p)>>(std::move(p)));
                            projectorMap_[decomposition_.id(*sd)][decomposition_.id(*mortar)] = projectors_.size() - 1;
                        } else {
                            DUNE_THROW(Dune::InvalidStateException, "Unexpected trace grid visited");
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

    // TODO: store mapping (use map with pointers as hashes?) to avoid many iterations
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
};

/*!
 * \ingroup MultiDomain
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

    //! Create a model from all inserted mortars & subdomains using default projectors
    Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...> make() const
    {
        return {
            decompositionFactory_.make(),
            solvers_,
            [] (const auto& mortarGG, const auto&, const auto& traceGrid, const auto&) {
                return FVDefaultProjector<MortarSolutionVector>{mortarGG, traceGrid};
            }
        };
    }

    //! Create a model from all inserted mortars & subdomains with a custom projector factory
    template<std::invocable T>
    Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...> make(T&& t) const
    { return {decompositionFactory_.make(), solvers_, std::forward<T>(t)}; }

 private:
    DecompositionFactory<MortarGridGeometry, SubDomainGridGeometries...> decompositionFactory_;
    std::vector<SolverVariant> solvers_;
};

} // end namespace Dumux::Mortar

#endif
