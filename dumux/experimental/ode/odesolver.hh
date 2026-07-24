// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Infrastructure for solving ordinary differential equations with multi-stage methods.
 */
#ifndef DUMUX_TIMESTEPPING_ODE_SOLVER_HH
#define DUMUX_TIMESTEPPING_ODE_SOLVER_HH

#include <algorithm>
#include <cmath>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/typetraits.hh>

#include <dumux/assembly/partialreassembler.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/variablesbackend.hh>
#include <dumux/experimental/common/variables.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux::Experimental::Detail::ODE {

template<class ODESystem>
struct VariablesChooser
{ using Type = Dumux::Experimental::Variables<typename ODESystem::SolutionVector>; };

template<class ODESystem>
    requires requires { typename ODESystem::Variables; }
struct VariablesChooser<ODESystem>
{ using Type = typename ODESystem::Variables; };

template<class ODESystem>
using Variables = typename VariablesChooser<ODESystem>::Type;

template<class To, class From>
void assign(To& to, const From& from)
{
    if constexpr (std::is_assignable_v<To&, From>)
        to = from;
    else
        static_assert(AlwaysFalse<To>::value, "Cannot assign ODE quantity to the requested output type.");
}

template<class Value>
void setZero(Value& value)
{
    value = 0.0;
}

template<class Matrix>
void setIdentity(Matrix& matrix)
{
    if constexpr (Dune::IsNumber<Matrix>::value)
        matrix = 1.0;
    else
    {
        matrix = 0.0;
        for (std::size_t i = 0; i < matrix.N(); ++i)
            matrix[i][i] = 1.0;
    }
}

template<class Vector, class Scalar>
void addScaled(const Scalar factor, const Vector& x, Vector& y)
{
    Dumux::DofBackend<Vector>::axpy(factor, x, y);
}

template<class Matrix, class Scalar>
void scaleMatrix(Matrix& matrix, const Scalar factor)
{
    matrix *= factor;
}

template<class Matrix, class Scalar>
void addScaledMatrix(const Scalar factor, const Matrix& x, Matrix& y)
{
    if constexpr (Dune::IsNumber<Matrix>::value)
        y += factor*x;
    else
    {
        for (std::size_t row = 0; row < y.N(); ++row)
            for (std::size_t col = 0; col < y.M(); ++col)
                y[row][col] += factor*x[row][col];
    }
}

template<class ODESystem, class Vars, class Residual>
void evaluateRhs(const ODESystem& odeSystem, const Vars& vars, Residual& rhs)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Residual& r) { system.rhs(v, r); })
        odeSystem.rhs(vars, rhs);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.rhs(v); })
        assign(rhs, odeSystem.rhs(vars));
    else
        static_assert(AlwaysFalse<ODESystem>::value,
            "ODE systems must provide rhs(vars, rhs) or rhs(vars).");
}

template<class ODESystem, class Vars, class Residual>
void evaluateStorage(const ODESystem& odeSystem, const Vars& vars, Residual& storage)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Residual& s) { system.storage(v, s); })
        odeSystem.storage(vars, storage);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.storage(v); })
        assign(storage, odeSystem.storage(vars));
    else
        assign(storage, Dumux::VariablesBackend<Vars>::dofs(vars));
}

template<class ODESystem, class Vars, class Jacobian>
void evaluateRhsJacobian(const ODESystem& odeSystem, const Vars& vars, Jacobian& jacobian)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Jacobian& j) { system.rhsJacobian(v, j); })
        odeSystem.rhsJacobian(vars, jacobian);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.rhsJacobian(v); })
        assign(jacobian, odeSystem.rhsJacobian(vars));
    else if constexpr (requires(const ODESystem& system, const Vars& v, Jacobian& j) { system.jacobian(v, j); })
        odeSystem.jacobian(vars, jacobian);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.jacobian(v); })
        assign(jacobian, odeSystem.jacobian(vars));
    else
        DUNE_THROW(Dune::NotImplemented,
            "Implicit ODE stages require rhsJacobian(vars, jacobian), rhsJacobian(vars), jacobian(vars, jacobian), or jacobian(vars).");
}

template<class ODESystem, class Vars, class Residual, class Jacobian>
void evaluateStorageJacobian(const ODESystem& odeSystem, const Vars& vars, Jacobian& jacobian)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Jacobian& j) { system.storageJacobian(v, j); })
        odeSystem.storageJacobian(vars, jacobian);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.storageJacobian(v); })
        assign(jacobian, odeSystem.storageJacobian(vars));
    else if constexpr (requires(const ODESystem& system, const Vars& v, Residual& s) { system.storage(v, s); }
                       || requires(const ODESystem& system, const Vars& v) { system.storage(v); })
        DUNE_THROW(Dune::NotImplemented,
            "ODE systems with custom storage require storageJacobian(vars, jacobian) or storageJacobian(vars).");
    else
        setIdentity(jacobian);
}

template<class Vars, class TimeLevel>
void updateTime(Vars& vars, const TimeLevel& timeLevel)
{
    if constexpr (requires(Vars& v, const TimeLevel& t) { v.updateTime(t); })
        vars.updateTime(timeLevel);
    else
        static_assert(AlwaysFalse<Vars>::value,
            "ODE variables must provide updateTime(timeLevel) or use Dumux::Experimental::Variables.");
}

template<class Vector>
auto norm(const Vector& vector)
{
    if constexpr (Dune::IsNumber<Vector>::value)
    {
        using std::abs;
        return abs(vector);
    }
    else if constexpr (requires(const Vector& v) { v.two_norm(); })
        return vector.two_norm();
    else
    {
        using Scalar = std::decay_t<decltype(norm(*std::begin(vector)))>;
        Scalar result = 0.0;
        for (const auto& entry : vector)
        {
            const auto entryNorm = norm(entry);
            result += entryNorm*entryNorm;
        }

        using std::sqrt;
        return sqrt(result);
    }
}

template<class Matrix, class Vector>
bool solve(const Matrix& matrix, Vector& x, const Vector& b)
{
    if constexpr (Dune::IsNumber<Matrix>::value)
    {
        x = b/matrix;
        return true;
    }
    else if constexpr (requires(const Matrix& A, Vector& y, const Vector& rhs) { A.solve(y, rhs); })
    {
        matrix.solve(x, b);
        return true;
    }
    else if constexpr (requires(Matrix A, const Vector& rhs, Vector& y) { A.invert(); A.mv(rhs, y); })
    {
        auto inverse = matrix;
        inverse.invert();
        inverse.mv(b, x);
        return true;
    }
    else
        static_assert(AlwaysFalse<Matrix>::value,
            "The default ODE linear solver only supports scalars and local dense matrices.");
}

} // end namespace Dumux::Experimental::Detail::ODE

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \brief Minimal linear solver for scalar and local dense ODE stage systems.
 *
 * This is intended for ordinary differential equations with scalar or small
 * fixed-size dense unknowns. Larger systems can pass any linear solver exposing
 * the usual DuMux Newton interface.
 */
template<class JacobianMatrix, class ResidualVector>
class ODELocalLinearSolver
{
public:
    void setResidualReduction(double) {}

    bool solve(const JacobianMatrix& matrix,
               ResidualVector& x,
               const ResidualVector& b) const
    { return Detail::ODE::solve(matrix, x, b); }

    auto norm(const ResidualVector& residual) const
    { return Detail::ODE::norm(residual); }
};

/*!
 * \ingroup Experimental
 * \brief Assembler that maps an ODE right-hand-side interface to the DuMux Newton interface.
 *
 * ODE systems are expected in the form \f$\dot u = f(u,t)\f$ and must export
 * `Scalar`, `SolutionVector`, `ResidualType`, and `JacobianMatrix`. They provide
 * `rhs(vars, rhs)` or `rhs(vars)`. Implicit stages additionally require
 * `jacobian(vars, jacobian)`/`jacobian(vars)` or the equivalent
 * `rhsJacobian` overloads. Optionally, systems can provide
 * `storage(vars, storage)` and `storageJacobian(vars, jacobian)` for equations
 * \f$\partial_t M(u,t) = f(u,t)\f$.
 */
template<class ODESystem>
class MultiStageODEAssembler
{
public:
    using Scalar = typename ODESystem::Scalar;
    using SolutionVector = typename ODESystem::SolutionVector;
    using ResidualType = typename ODESystem::ResidualType;
    using JacobianMatrix = typename ODESystem::JacobianMatrix;
    using Variables = Detail::ODE::Variables<ODESystem>;
    using StageParams = MultiStageParams<Scalar>;

    explicit MultiStageODEAssembler(std::shared_ptr<const ODESystem> odeSystem)
    : odeSystem_(std::move(odeSystem))
    {}

    explicit MultiStageODEAssembler(const ODESystem& odeSystem)
    : MultiStageODEAssembler(std::make_shared<ODESystem>(odeSystem))
    {}

    void setLinearSystem() {}

    JacobianMatrix& jacobian()
    { return jacobian_; }

    const JacobianMatrix& jacobian() const
    { return jacobian_; }

    ResidualType& residual()
    { return residual_; }

    const ResidualType& residual() const
    { return residual_; }

    void assembleResidual(const Variables& vars)
    {
        if (!stageParams_)
            DUNE_THROW(Dune::InvalidStateException, "No ODE stage prepared before residual assembly");

        if (stageParams_->size() != rhsTerms_.size() || stageParams_->size() != storageTerms_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of ODE stage residuals");

        evaluateTerms_(vars, storageTerms_.back(), rhsTerms_.back());

        residual_ = zeroResidualLike_(vars);
        for (std::size_t k = 0; k < stageParams_->size(); ++k)
        {
            if (!stageParams_->skipTemporal(k))
                Detail::ODE::addScaled(stageParams_->temporalWeight(k), storageTerms_[k], residual_);
            if (!stageParams_->skipSpatial(k))
                Detail::ODE::addScaled(-stageParams_->spatialWeight(k), rhsTerms_[k], residual_);
        }
    }

    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);

        const auto curStage = stageParams_->size() - 1;
        Detail::ODE::evaluateStorageJacobian<ODESystem, Variables, ResidualType>(*odeSystem_, vars, jacobian_);
        Detail::ODE::scaleMatrix(jacobian_, stageParams_->temporalWeight(curStage));

        if (!stageParams_->skipSpatial(curStage))
        {
            auto rhsJacobian = jacobian_;
            Detail::ODE::setZero(rhsJacobian);
            Detail::ODE::evaluateRhsJacobian(*odeSystem_, vars, rhsJacobian);
            Detail::ODE::addScaledMatrix(-stageParams_->spatialWeight(curStage), rhsJacobian, jacobian_);
        }
    }

    void prepareStage(Variables& vars,
                      std::shared_ptr<const StageParams> params)
    {
        const auto curStage = params->size() - 1;
        const auto prevStage = stageParams_ ? stageParams_->size() - 1 : 0;
        if (curStage != prevStage + 1)
            DUNE_THROW(Dune::InvalidStateException,
                "Can only prepare ODE stages consecutively (current stage: " << curStage
                << ", previous stage: " << prevStage << ")");

        stageParams_ = std::move(params);

        if (curStage == 1)
        {
            previousSolution_ = Dumux::VariablesBackend<Variables>::dofs(vars);

            const auto t = stageParams_->timeAtStage(0);
            Detail::ODE::updateTime(vars, TimeLevel<Scalar>{t, t, stageParams_->timeStepFraction(0)});

            storageTerms_.push_back(zeroResidualLike_(vars));
            rhsTerms_.push_back(zeroResidualLike_(vars));
            evaluateTerms_(vars, storageTerms_.back(), rhsTerms_.back());
        }

        const auto t = stageParams_->timeAtStage(curStage);
        const auto prevT = stageParams_->timeAtStage(0);
        const auto dtFraction = stageParams_->timeStepFraction(curStage);
        Detail::ODE::updateTime(vars, TimeLevel<Scalar>{t, prevT, dtFraction});

        storageTerms_.push_back(zeroResidualLike_(vars));
        rhsTerms_.push_back(zeroResidualLike_(vars));
    }

    void clearStages()
    {
        storageTerms_.clear();
        rhsTerms_.clear();
        stageParams_.reset();
    }

    const SolutionVector& prevSol() const
    { return previousSolution_; }

    void resetTimeStep(const SolutionVector&)
    { clearStages(); }

    const ODESystem& odeSystem() const
    { return *odeSystem_; }

private:
    using Backend = Dumux::VariablesBackend<Variables>;
    using ResidualBackend = Dumux::DofBackend<ResidualType>;

    ResidualType zeroResidualLike_(const Variables& vars) const
    {
        auto result = ResidualBackend::zeros(Backend::size(Backend::dofs(vars)));
        Detail::ODE::setZero(result);
        return result;
    }

    void evaluateTerms_(const Variables& vars,
                        ResidualType& storage,
                        ResidualType& rhs) const
    {
        Detail::ODE::setZero(storage);
        Detail::ODE::setZero(rhs);
        Detail::ODE::evaluateStorage(*odeSystem_, vars, storage);
        Detail::ODE::evaluateRhs(*odeSystem_, vars, rhs);
    }

    std::shared_ptr<const ODESystem> odeSystem_;
    ResidualType residual_;
    JacobianMatrix jacobian_;
    SolutionVector previousSolution_;
    std::vector<ResidualType> storageTerms_;
    std::vector<ResidualType> rhsTerms_;
    std::shared_ptr<const StageParams> stageParams_;
};

/*!
 * \ingroup Experimental
 * \brief Convenience wrapper composing ODE assembly, Newton, and multi-stage stepping.
 */
template<class ODESystem,
         class LinearSolver = ODELocalLinearSolver<typename ODESystem::JacobianMatrix, typename ODESystem::ResidualType>,
         class Reassembler = DefaultPartialReassembler,
         class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator>>
class ODESolver
{
public:
    using Assembler = MultiStageODEAssembler<ODESystem>;
    using NonlinearSolver = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using Scalar = typename Assembler::Scalar;
    using SolutionVector = typename Assembler::SolutionVector;
    using Variables = typename Assembler::Variables;
    using Method = MultiStageMethod<Scalar>;
    using TimeStepper = MultiStageTimeStepper<NonlinearSolver, Scalar>;

    ODESolver(std::shared_ptr<const ODESystem> odeSystem,
              std::shared_ptr<const Method> method,
              const std::string& paramGroup = "")
    : ODESolver(std::move(odeSystem),
                std::make_shared<LinearSolver>(),
                std::move(method),
                Dune::MPIHelper::getCommunication(),
                paramGroup)
    {}

    ODESolver(std::shared_ptr<const ODESystem> odeSystem,
              std::shared_ptr<const Method> method,
              const Comm& comm,
              const std::string& paramGroup = "")
    : ODESolver(std::move(odeSystem),
                std::make_shared<LinearSolver>(),
                std::move(method),
                comm,
                paramGroup)
    {}

    ODESolver(std::shared_ptr<const ODESystem> odeSystem,
              std::shared_ptr<LinearSolver> linearSolver,
              std::shared_ptr<const Method> method,
              const std::string& paramGroup = "")
    : ODESolver(std::move(odeSystem),
                std::move(linearSolver),
                std::move(method),
                Dune::MPIHelper::getCommunication(),
                paramGroup)
    {}

    ODESolver(std::shared_ptr<const ODESystem> odeSystem,
              std::shared_ptr<LinearSolver> linearSolver,
              std::shared_ptr<const Method> method,
              const Comm& comm,
              const std::string& paramGroup = "")
    : odeSystem_(std::move(odeSystem))
    , assembler_(std::make_shared<Assembler>(odeSystem_))
    , linearSolver_(std::move(linearSolver))
    , nonlinearSolver_(std::make_shared<NonlinearSolver>(assembler_, linearSolver_, comm, paramGroup, "ODE", 0))
    , method_(std::move(method))
    , timeStepper_(nonlinearSolver_, method_, paramGroup)
    {}

    void step(Variables& vars, const Scalar t, const Scalar dt)
    { timeStepper_.step(vars, t, dt); }

    void step(Variables& vars, TimeLoopBase<Scalar>& timeLoop)
    { timeStepper_.step(vars, timeLoop); }

    void solve(Variables& vars, const Scalar endTime, const Scalar timeStepSize)
    {
        if (!(timeStepSize > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "ODE time step size has to be positive");

        auto time = vars.timeLevel().current();
        while (time < endTime)
        {
            using std::min;
            const auto dt = min(timeStepSize, endTime - time);
            step(vars, time, dt);
            time = vars.timeLevel().current();
        }
    }

    Assembler& assembler()
    { return *assembler_; }

    const Assembler& assembler() const
    { return *assembler_; }

    NonlinearSolver& nonlinearSolver()
    { return *nonlinearSolver_; }

    const NonlinearSolver& nonlinearSolver() const
    { return *nonlinearSolver_; }

    LinearSolver& linearSolver()
    { return *linearSolver_; }

    const LinearSolver& linearSolver() const
    { return *linearSolver_; }

private:
    std::shared_ptr<const ODESystem> odeSystem_;
    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<LinearSolver> linearSolver_;
    std::shared_ptr<NonlinearSolver> nonlinearSolver_;
    std::shared_ptr<const Method> method_;
    TimeStepper timeStepper_;
};

} // end namespace Dumux::Experimental

#endif
