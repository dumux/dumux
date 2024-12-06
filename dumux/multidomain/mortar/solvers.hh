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

#include <dumux/common/properties.hh>
#include <dumux/discretization/facetgrid.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "model.hh"
#include "solverinterface.hh"

namespace Dumux::Mortar {

//! Default solver for stationary single-domain problems, discretized by finite volume schemes
template<typename LinearSolver, typename TypeTag>
class DefaultSubDomainSolver : public SubDomainSolver<
    GetPropType<TypeTag, Properties::MortarSolutionVector>,
    GetPropType<TypeTag, Properties::MortarGrid>,
    GetPropType<TypeTag, Properties::GridGeometry>
>
{
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using MortarSolutionVector = GetPropType<TypeTag, Properties::MortarSolutionVector>;

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    using ParentType = SubDomainSolver<
        MortarSolutionVector,
        GetPropType<TypeTag, Properties::MortarGrid>,
        GetPropType<TypeTag, Properties::GridGeometry>
    >;

 public:
    using typename ParentType::GridGeometry;
    using typename ParentType::TraceGrid;

    DefaultSubDomainSolver(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType{std::move(gridGeometry)}
    {
        x_.resize(this->gridGeometry()->numDofs());
        x_ = 0.0;

        problem_ = std::make_shared<Problem>(this->gridGeometry());
        gridVariables_ = std::make_shared<GridVariables>(problem_, this->gridGeometry());
        gridVariables_->init(x_);

        assembler_ = std::make_shared<Assembler>(problem_, this->gridGeometry(), gridVariables_);
        newtonSolver_ = std::make_unique<NewtonSolver>(assembler_, std::make_shared<LinearSolver>());
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

} // end namespace Dumux::Mortar

#endif
