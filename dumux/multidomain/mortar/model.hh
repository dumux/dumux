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

    virtual ~Solver() = default;
    explicit Solver(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_{std::move(gridGeometry)}
    {}

    virtual void solve() = 0;
    // virtual void

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
        decomposition_.visitMortars([&] (const auto& gg) { numMortarDofs += gg->numDofs(); });
        mortarSolution_.resize(numMortarDofs);
    }

    Decomposition decomposition_;
    std::vector<SolverVariant> solvers_;
    MortarSolutionVector mortarSolution_;
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
