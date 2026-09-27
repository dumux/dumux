// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Default subdomain solver for single-domain problems.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_SOLVERS_HH
#define DUMUX_MULTIDOMAIN_MORTAR_SOLVERS_HH

#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

#include <dumux/common/properties.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "solverinterface.hh"
#include "couplingmanager.hh"
#include "properties.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief Default solver for stationary single-domain problems, discretized by finite volume
 *        schemes.
 *
 * The problem of the type tag is constructed as `Problem(gridGeometry, couplingManager,
 * paramGroup)`, receiving the SubDomainCouplingManager of the subdomain, which it asks where
 * mortar data is imposed and what it is. It returns the conjugate trace on the mortar domain
 * with a given id through `problem.assembleTraceVariables(mortarId, gridVariables, x)`,
 * usually by assembleTrace with an integrand of its physics. The type tag defines the
 * properties MortarGrid and MortarSolutionVector.
 *
 * The trace data is taken per trace vertex for control-volume finite element schemes and
 * per trace cell otherwise, as the coupling manager declares; setTraceDataOrder changes
 * that before the model is built.
 */
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

    using GridView = typename ParentType::GridGeometry::GridView;
    using Communication = std::decay_t<decltype(std::declval<const GridView&>().comm())>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using SolverTraits = Dumux::LinearSolverTraits<typename ParentType::GridGeometry>;
    using LinearSolver = IstlSolverFactoryBackend<SolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver, PartialReassembler<Assembler>, Communication>;

 public:
    using typename ParentType::GridGeometry;
    using typename ParentType::Trace;
    using typename ParentType::MortarSolutionVector;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using MortarCouplingManager = SubDomainCouplingManager<GridGeometry, GetPropType<TypeTag, Properties::MortarGrid>, MortarSolutionVector>;

    /*!
     * \brief A solver building the problem, the grid variables, the assembler and the Newton
     *        solver of the subdomain with the given grid geometry.
     * \param gridGeometry The grid geometry of the subdomain
     * \param paramGroup The parameter group of the problem and the solvers
     */
    DefaultSubDomainSolver(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType{std::move(gridGeometry)}
    {
        x_.resize(this->gridGeometry()->numDofs());
        x_ = 0.0;

        mortarCouplingManager_ = std::make_shared<MortarCouplingManager>(this->gridGeometry());
        problem_ = std::make_shared<Problem>(this->gridGeometry(), mortarCouplingManager_, paramGroup);
        gridVariables_ = std::make_shared<GridVariables>(problem_, this->gridGeometry());
        gridVariables_->init(x_);

        assembler_ = std::make_shared<Assembler>(problem_, this->gridGeometry(), gridVariables_);
        newtonSolver_ = std::make_unique<NewtonSolver>(assembler_, makeLinearSolver_(paramGroup),
                                                       this->gridGeometry()->gridView().comm(), paramGroup);
    }

    void solve() override
    { newtonSolver_->solve(x_); }

    void setTraceVariables(std::size_t mortarId, MortarSolutionVector trace) override
    { mortarCouplingManager_->setTraceVariables(mortarId, std::move(trace)); }

    void registerMortarTrace(std::shared_ptr<const Trace> trace, std::size_t mortarId) override
    { mortarCouplingManager_->registerTrace(trace, mortarId); }

    MortarSolutionVector assembleTraceVariables(std::size_t mortarId) const override
    { return problem_->assembleTraceVariables(mortarId, *gridVariables_, x_); }

    void setCouplingMode(CouplingMode mode) override
    { mortarCouplingManager_->setCouplingMode(mode); }

    void setHomogeneous(bool homogeneous) override
    { mortarCouplingManager_->setHomogeneous(homogeneous); }

    bool isFloating() const override
    { return mortarCouplingManager_->isFloating(); }

    std::size_t traceDataOrder() const override
    { return mortarCouplingManager_->traceDataOrder(); }

    //! Set the order of the trace data this subdomain accepts, before the model is built
    void setTraceDataOrder(std::size_t order)
    { mortarCouplingManager_->setTraceDataOrder(order); }

    //! The subdomain problem
    Problem& problem() { return *problem_; }
    //! The subdomain problem
    const Problem& problem() const { return *problem_; }

    //! The grid variables of the subdomain
    GridVariables& gridVariables() { return *gridVariables_; }
    //! The grid variables of the subdomain
    const GridVariables& gridVariables() const { return *gridVariables_; }

    //! The solution of the last subdomain solve
    SolutionVector& solution() { return x_; }
    //! The solution of the last subdomain solve
    const SolutionVector& solution() const { return x_; }

 private:
    //! the linear solver takes its communicator from the grid view where the grid can communicate
    std::shared_ptr<LinearSolver> makeLinearSolver_(const std::string& paramGroup) const
    {
        if constexpr (SolverTraits::canCommunicate)
            return std::make_shared<LinearSolver>(this->gridGeometry()->gridView(), this->gridGeometry()->dofMapper(), paramGroup);
        else
            return std::make_shared<LinearSolver>(paramGroup);
    }

    std::shared_ptr<MortarCouplingManager> mortarCouplingManager_;
    std::shared_ptr<GridVariables> gridVariables_;
    std::shared_ptr<Problem> problem_;
    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<NewtonSolver> newtonSolver_;
    SolutionVector x_;
};

} // end namespace Dumux::Mortar

#endif
