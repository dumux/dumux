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
 * \brief Default subdomain solver implementations.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_SOLVERS_HH
#define DUMUX_MULTIDOMAIN_MORTAR_SOLVERS_HH

#include <memory>
#include <utility>
#include <tuple>

#include <dumux/linear/istlsolverfactorybackend.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/facetgrid.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "solverinterface.hh"
#include "properties.hh"

namespace Dumux::Mortar {

//! Default solver for stationary single-domain problems, discretized by finite volume schemes
template<typename TypeTag>
class DefaultSubDomainSolver : public SubDomainSolver<
    GetPropType<TypeTag, Properties::MortarSolutionVector>,
    GetPropType<TypeTag, Properties::MortarGrid>,
    GetPropType<TypeTag, Properties::GridGeometry>
>
{
    using ParentType = SubDomainSolver<
        GetPropType<TypeTag, Properties::MortarSolutionVector>,
        GetPropType<TypeTag, Properties::MortarGrid>,
        GetPropType<TypeTag, Properties::GridGeometry>
    >;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = IstlSolverFactoryBackend<
        Dumux::LinearSolverTraits<GetPropType<TypeTag, Properties::GridGeometry>>,
        LinearAlgebraTraitsFromAssembler<Assembler>
    >;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

 public:
    using typename ParentType::GridGeometry;
    using typename ParentType::TraceGrid;
    using typename ParentType::MortarSolutionVector;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    DefaultSubDomainSolver(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType{std::move(gridGeometry)}
    {
        x_.resize(this->gridGeometry()->numDofs());
        x_ = 0.0;

        problem_ = std::make_shared<Problem>(this->gridGeometry(), paramGroup);
        gridVariables_ = std::make_shared<GridVariables>(problem_, this->gridGeometry());
        gridVariables_->init(x_);

        assembler_ = std::make_shared<Assembler>(problem_, this->gridGeometry(), gridVariables_);
        newtonSolver_ = std::make_unique<NewtonSolver>(assembler_, std::make_shared<LinearSolver>(paramGroup));
    }

    void solve() override
    { newtonSolver_->solve(x_); }

    void setTraceVariables(std::size_t mortarId, MortarSolutionVector trace) override
    { problem_->setTraceVariables(mortarId, std::move(trace)); }

    void registerMortarTrace(std::shared_ptr<const TraceGrid> traceGrid, std::size_t mortarId) override
    { problem_->registerMortarTrace(traceGrid, mortarId); }

    MortarSolutionVector assembleTraceVariables(std::size_t mortarId) const override
    { return problem_->assembleTraceVariables(mortarId, *gridVariables_, x_); }

    Problem& problem() { return *problem_; }
    const Problem& problem() const { return *problem_; }

    GridVariables& gridVariables() { return *gridVariables_; }
    const GridVariables& gridVariables() const { return *gridVariables_; }

    SolutionVector& solution() { return x_; }
    const SolutionVector& solution() const { return x_; }

 private:
    std::shared_ptr<GridVariables> gridVariables_;
    std::shared_ptr<Problem> problem_;
    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<NewtonSolver> newtonSolver_;
    SolutionVector x_;
};

//! Default solver for stationary single-domain problems, discretized by the staggered grid scheme
template<typename MomentumTypeTag, typename MassTypeTag>
class DefaultStaggeredFreeFlowSubDomainSolver : public SubDomainSolver<
    GetPropType<MomentumTypeTag, Properties::MortarSolutionVector>,
    GetPropType<MomentumTypeTag, Properties::MortarGrid>,
    GetPropType<MomentumTypeTag, Properties::GridGeometry>
>
{
    using MortarSolutionVector = GetPropType<MomentumTypeTag, Properties::MortarSolutionVector>;
    using ParentType = SubDomainSolver<
        MortarSolutionVector,
        GetPropType<MomentumTypeTag, Properties::MortarGrid>,
        GetPropType<MomentumTypeTag, Properties::GridGeometry>
    >;

    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    using MDTraits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename MDTraits::SolutionVector;

 public:
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;

    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;

    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;

 private:
    using Assembler = MultiDomainFVAssembler<MDTraits, CouplingManager, DiffMethod::numeric>;
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;

    static constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    static constexpr auto massIdx = CouplingManager::freeFlowMassIndex;

 public:
    using typename ParentType::GridGeometry;
    using typename ParentType::TraceGrid;

    DefaultStaggeredFreeFlowSubDomainSolver(std::shared_ptr<const MomentumGridGeometry> momentumGridGeometry,
                                            std::shared_ptr<const MassGridGeometry> massGridGeometry)
    : ParentType{momentumGridGeometry}
    , momentumGridGeometry_{momentumGridGeometry}
    , massGridGeometry_{massGridGeometry}
    {
        couplingManager_ = std::make_shared<CouplingManager>();
        momentumProblem_ = std::make_shared<MomentumProblem>(momentumGridGeometry_, couplingManager_);
        massProblem_ = std::make_shared<MassProblem>(massGridGeometry_, couplingManager_);

        x_[momentumIdx].resize(momentumGridGeometry->numDofs());
        x_[massIdx].resize(massGridGeometry->numDofs());
        x_ = 0.0;

        momentumGridVariables_ = std::make_shared<MomentumGridVariables>(momentumProblem_, momentumGridGeometry_);
        massGridVariables_ = std::make_shared<MassGridVariables>(massProblem_, massGridGeometry_);
        couplingManager_->init(momentumProblem_, massProblem_, std::make_tuple(momentumGridVariables_, massGridVariables_), x_);
        massGridVariables_->init(x_[massIdx]);
        momentumGridVariables_->init(x_[momentumIdx]);

        assembler_ = std::make_shared<Assembler>(std::make_tuple(momentumProblem_, massProblem_),
                                                 std::make_tuple(momentumGridGeometry_, massGridGeometry_),
                                                 std::make_tuple(momentumGridVariables_, massGridVariables_),
                                                 couplingManager_);
        nonLinearSolver_ = std::make_shared<NewtonSolver>(assembler_, std::make_shared<LinearSolver>(), couplingManager_);
    }

    void solve() override
    { nonLinearSolver_->solve(x_); }

    void setTraceVariables(std::size_t mortarId, MortarSolutionVector trace) override
    { momentumProblem_->setTraceVariables(mortarId, std::move(trace)); }

    void registerMortarTrace(std::shared_ptr<const TraceGrid> traceGrid, std::size_t mortarId) override
    { momentumProblem_->registerMortarTrace(traceGrid, mortarId); }

    MortarSolutionVector assembleTraceVariables(std::size_t mortarId) const override
    { return momentumProblem_->assembleTraceVariables(mortarId, *momentumGridVariables_, x_[momentumIdx]); }

    MomentumProblem& momentumProblem() { return *momentumProblem_; }
    const MomentumProblem& momentumProblem() const { return *momentumProblem_; }

    MassProblem& massProblem() { return *massProblem_; }
    const MassProblem& massProblem() const { return *massProblem_; }

    MomentumGridVariables& momentumGridVariables() { return *momentumGridVariables_; }
    const MomentumGridVariables& momentumGridVariables() const { return *momentumGridVariables_; }

    MassGridVariables& massGridVariables() { return *massGridVariables_; }
    const MassGridVariables& massGridVariables() const { return *massGridVariables_; }

    SolutionVector& solution() { return x_; }
    const SolutionVector& solution() const { return x_; }

    const auto& massSolution() const { return x_[massIdx]; }
    const auto& momentumSolution() const { return x_[momentumIdx]; }

 private:
    std::shared_ptr<const MomentumGridGeometry> momentumGridGeometry_;
    std::shared_ptr<const MassGridGeometry> massGridGeometry_;

    std::shared_ptr<MomentumProblem> momentumProblem_;
    std::shared_ptr<MassProblem> massProblem_;

    std::shared_ptr<MomentumGridVariables> momentumGridVariables_;
    std::shared_ptr<MassGridVariables> massGridVariables_;

    std::shared_ptr<CouplingManager> couplingManager_;
    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<LinearSolver> linearSolver_;
    std::shared_ptr<NewtonSolver> nonLinearSolver_;

    SolutionVector x_;
};

} // end namespace Dumux::Mortar

#endif
