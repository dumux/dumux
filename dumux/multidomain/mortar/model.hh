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

#include "decomposition.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Abstract base class for subdomain solvers.
 */
template<typename MortarSolutionVector, typename GG>
class Solver
{
 public:
    using GridGeometry = GG;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

    virtual ~Solver() = default;
    explicit Solver(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_{std::move(gridGeometry)}
    {}

    //! Solve the subdomain problem
    virtual void solve() = 0;

    //! Set the mortar boundary condition for the mortar with the given id
    virtual void setMortar(std::size_t, MortarSolutionVector) = 0;

    //! Register a mapping from the given sub-control volume face to the mortar domain with the given id
    virtual void registerMortarScvf(const Element&, const SubControlVolumeFace&, std::size_t) = 0;

    //! Assemble the variables on the trace overlapping with the given mortar domain
    virtual MortarSolutionVector assembleTraceVariables(std::size_t) const = 0;

    //! Return the underlying grid geometry
    const std::shared_ptr<const GridGeometry>& gridGeometry() const
    { return gridGeometry_; }

 private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
};

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
 public:
    using Decomposition = Mortar::Decomposition<MortarGridGeometry, SubDomainGridGeometries...>;
    using SolverVariant = std::variant<std::shared_ptr<Solver<MortarSolutionVector, SubDomainGridGeometries>>...>;

    void setMortar(const MortarSolutionVector& x)
    {
        decomposition_->visitMortar([&] (const auto& mortar) {
            const auto mortarId = decomposition_.id(*mortar);
            decomposition_->visitCoupledSubDomainsOf(*mortar, [&] (const auto& subDomain) {
                visitSolverFor_(subDomain, [&] (const auto& solver) {
                    MortarSolutionVector restricted(mortar->numDofs());
                    for (std::size_t i = 0; i < mortar->numDofs(); ++i)
                        restricted[i] = x[mortarDofOffsets_[mortarId] + i];
                    solver.setMortar(mortarId, restricted);
                });
            });
        });
    }

    void solveSubDomains()
    {
        // TODO: parallel
        std::for_each(solvers_.begin(), solvers_.end(), [] (auto& solver) {
            solver.solve();
        });
    }

    void assembleMortarResidual(MortarSolutionVector& residual) const
    {
        residual.resize(mortarSolution_.size());
        decomposition_->visitMortar([&] (const auto& mortar) {
            const auto mortarId = decomposition_.id(*mortar);
            decomposition_->visitCoupledSubDomainsOf(*mortar, [&] (const auto& subDomain) {
                visitSolverFor_(subDomain, [&] (const auto& solver) {
                    [[maybe_unused]] const auto vars = solver.assembleTraceVariables(mortarId);
                    // TODO: project
                    // if (projected.size() != mortar.numDofs())
                    //     DUNE_THROW(Dune::InvalidStateException, "Trace does not have the expected number of entries.");
                    // std::ranges::for_each(projected, [i=std::size_t{0}] (const auto& entry) mutable {
                    //     residual[mortarDofOffsets_[mortarId] + i++] += entry;
                    // });
                });
            });
        });
    }

    const MortarSolutionVector& mortarSolution() const
    { return mortarSolution_; }

 private:
    friend ModelFactory<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>;

    Model(Decomposition&& decomposition, std::vector<SolverVariant> solvers)
    : decomposition_{std::move(decomposition)}
    , solvers_{std::move(solvers)}
    {
        if (decomposition_.numberOfSubDomains() != solvers_.size())
            DUNE_THROW(Dune::InvalidStateException, "Number of solvers and subdomains do not match.");

        std::size_t numMortarDofs = 0;
        mortarDofOffsets_.resize(decomposition_.numberOfMortars());
        decomposition_.visitMortars([&] (const auto& mortar) {
            mortarDofOffsets_[decomposition_.id(*mortar)] = mortar->numDofs();
            numMortarDofs += mortar->numDofs();
        });
        std::exclusive_scan(
            mortarDofOffsets_.begin(), mortarDofOffsets_.end(),
            mortarDofOffsets_.begin(), std::size_t{0}
        );
        mortarSolution_.resize(numMortarDofs);
        std::ranges::fill(mortarSolution_, 0);
    }

    // TODO: store mapping (use map with pointers as hashes?) to avoid many iterations
    template<typename GridGeometry, typename Visitor>
    void visitSolverFor_(const GridGeometry& subDomain, Visitor&& v) const
    {
        for (const auto& variant : solvers_)
            if(
                std::visit([&] (const auto& solver) {
                    if (solver.gridGeometry().get() == &subDomain) { v(solver); return true; }
                    return false;
                }, variant)
            )
                return;
        DUNE_THROW(Dune::InvalidStateException, "Could not find solver matching the given subdomain");
    }

    Decomposition decomposition_;
    std::vector<SolverVariant> solvers_;
    MortarSolutionVector mortarSolution_;
    std::vector<std::size_t> mortarDofOffsets_;
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
    using SolverVariant = typename Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...>::SolverVariant;

    template<typename T>
    static constexpr bool isSupportedSolver
        = std::disjunction_v<std::is_same<typename T::GridGeometry, SubDomainGridGeometries>...>
        and std::derived_from<T, Solver<MortarSolutionVector, typename T::GridGeometry>>;

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

    //! Create a model from all inserted mortars & subdomains
    Model<MortarSolutionVector, MortarGridGeometry, SubDomainGridGeometries...> make() const
    { return {decompositionFactory_.make(), solvers_}; }

 private:
    DecompositionFactory<MortarGridGeometry, SubDomainGridGeometries...> decompositionFactory_;
    std::vector<SolverVariant> solvers_;
};

} // end namespace Dumux::Mortar

#endif
