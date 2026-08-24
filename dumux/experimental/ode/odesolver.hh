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
#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/variablesbackend.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux::Experimental::Detail::ODE {

template<class T, bool indexable = Dune::IsIndexable<T>::value>
struct ScalarType
{ using Type = T; };

template<class T>
struct ScalarType<T, true>
{
    using Element = std::decay_t<decltype(std::declval<T>()[0])>;
    using Type = typename ScalarType<Element>::Type;
};

template<class X>
using Scalar = typename ScalarType<X>::Type;

} // end namespace Dumux::Experimental::Detail::ODE

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \brief Represents a level of the independent variable during ODE integration
 *
 * \tparam Scalar the scalar type
 */
template<class Scalar>
class IndependentVariableLevel
{
public:
    /*!
     * \brief Constructs a level at a given value of the independent variable
     *
     * \param current  the current value of the independent variable
     */
    explicit IndependentVariableLevel(const Scalar current)
    : current_(current)
    , previous_(current)
    , stepFraction_(1.0)
    {}

    /*!
     * \brief Constructs a level occurring within an integration step
     *
     * \param current       the current value of the independent variable
     * \param previous      the value at the beginning of the integration step
     * \param stepFraction  the fraction of the integration step represented by this level
     */
    IndependentVariableLevel(const Scalar current,
                             const Scalar previous,
                             const Scalar stepFraction)
    : current_(current)
    , previous_(previous)
    , stepFraction_(stepFraction)
    {}

    /*!
     * \brief Returns the current value of the independent variable
     *
     * \return the current value of the independent variable
     */
    Scalar current() const
    { return current_; }

    /*!
     * \brief Returns the value at the beginning of the integration step
     *
     * \return the previous value of the independent variable
     */
    Scalar previous() const
    { return previous_; }

    /*!
     * \brief Returns the fraction of the integration step represented by this level
     *
     * \return the integration-step fraction
     */
    Scalar stepFraction() const
    { return stepFraction_; }

private:
    Scalar current_;
    Scalar previous_;
    Scalar stepFraction_;
};

/*!
 * \ingroup Experimental
 * \brief Stores an ODE solution and its independent-variable level. Generalization of `experimental/common/variables` which are tailored to time as the independent variable.
 *
 * \tparam X the solution-vector type
 */
template<class X>
class ODEVariables
{
public:
    using SolutionVector = X;
    using Scalar = Detail::ODE::Scalar<X>;
    using Level = IndependentVariableLevel<Scalar>;

    /*!
     * \brief Constructs default-initialized ODE variables
     */
    ODEVariables()
    : solution_()
    , level_(0.0)
    {}

    /*!
     * \brief Constructs ODE variables from a solution
     *
     * \param solution initial solution
     * \param level    initial independent-variable level
     */
    explicit ODEVariables(const SolutionVector& solution,
                          const Level& level = Level{0.0})
    : solution_(solution)
    , level_(level)
    {}

    /*!
     * \brief Constructs ODE variables by moving a solution
     *
     * \param solution initial solution
     * \param level    initial independent-variable level
     */
    explicit ODEVariables(SolutionVector&& solution,
                          const Level& level = Level{0.0})
    : solution_(std::move(solution))
    , level_(level)
    {}

    /*!
     * \brief Returns the independent-variable level
     *
     * \return the independent-variable level
     */
    const Level& independentVariableLevel() const
    { return level_; }

    /*!
     * \brief Returns the solution degrees of freedom
     *
     * \return the solution degrees of freedom
     */
    const SolutionVector& dofs() const
    { return solution_; }

    /*!
     * \brief Returns the solution degrees of freedom
     *
     * \return the solution degrees of freedom
     */
    SolutionVector& dofs()
    { return solution_; }

    /*!
     * \brief Updates the solution
     *
     * \param solution the new solution
     */
    void update(const SolutionVector& solution)
    { solution_ = solution; }

    /*!
     * \brief Updates the independent-variable level
     *
     * \param level the new independent-variable level
     */
    void updateIndependentVariable(const Level& level)
    { level_ = level; }

    /*!
     * \brief Updates the solution and independent-variable level
     *
     * \param solution new solution
     * \param level    new independent-variable level
     */
    void update(const SolutionVector& solution, const Level& level)
    {
        solution_ = solution;
        level_ = level;
    }

private:
    SolutionVector solution_;
    Level level_;
};

} // end namespace Dumux::Experimental

namespace Dumux::Experimental::Detail::ODE {

template<class ODESystem>
struct VariablesChooser
{ using Type = Dumux::Experimental::ODEVariables<typename ODESystem::SolutionVector>; };

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
void evaluateDerivativeQuantity(const ODESystem& odeSystem, const Vars& vars, Residual& value)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Residual& result) { system.derivativeQuantity(v, result); })
        odeSystem.derivativeQuantity(vars, value);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.derivativeQuantity(v); })
        assign(value, odeSystem.derivativeQuantity(vars));
    else
        assign(value, Dumux::VariablesBackend<Vars>::dofs(vars));
}

template<class ODESystem, class Vars, class Jacobian>
void evaluateRhsJacobian(const ODESystem& odeSystem, const Vars& vars, Jacobian& jacobian)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Jacobian& j) { system.rhsJacobian(v, j); })
        odeSystem.rhsJacobian(vars, jacobian);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.rhsJacobian(v); })
        assign(jacobian, odeSystem.rhsJacobian(vars));
    else
        DUNE_THROW(Dune::NotImplemented,
            "Implicit ODE stages require rhsJacobian(vars, jacobian) or rhsJacobian(vars).");
}

template<class ODESystem, class Vars, class Residual, class Jacobian>
void evaluateDerivativeQuantityJacobian(const ODESystem& odeSystem, const Vars& vars, Jacobian& jacobian)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Jacobian& j) { system.derivativeQuantityJacobian(v, j); })
        odeSystem.derivativeQuantityJacobian(vars, jacobian);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.derivativeQuantityJacobian(v); })
        assign(jacobian, odeSystem.derivativeQuantityJacobian(vars));
    else if constexpr (requires(const ODESystem& system, const Vars& v, Residual& result) { system.derivativeQuantity(v, result); }
                       || requires(const ODESystem& system, const Vars& v) { system.derivativeQuantity(v); })
        DUNE_THROW(Dune::NotImplemented,
            "ODE systems with a custom derivative quantity require derivativeQuantityJacobian(vars, jacobian) or derivativeQuantityJacobian(vars).");
    else
        setIdentity(jacobian);
}

template<class Vars, class Level>
void updateIndependentVariable(Vars& vars, const Level& level)
{
    if constexpr (requires(Vars& v, const Level& l) { v.updateIndependentVariable(l); })
        vars.updateIndependentVariable(level);
    else
        static_assert(AlwaysFalse<Vars>::value,
            "ODE variables must provide updateIndependentVariable(level) or use Dumux::Experimental::ODEVariables.");
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

        if (stageParams_->size() != rhsTerms_.size() || stageParams_->size() != derivativeQuantities_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of ODE stage residuals");

        evaluateTerms_(vars, derivativeQuantities_.back(), rhsTerms_.back());

        residual_ = zeroResidualLike_(vars);
        for (std::size_t k = 0; k < stageParams_->size(); ++k)
        {
            if (!stageParams_->skipTemporal(k))
                Detail::ODE::addScaled(stageParams_->temporalWeight(k), derivativeQuantities_[k], residual_);
            if (!stageParams_->skipSpatial(k))
                Detail::ODE::addScaled(-stageParams_->spatialWeight(k), rhsTerms_[k], residual_);
        }
    }

    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);

        const auto curStage = stageParams_->size() - 1;
        Detail::ODE::evaluateDerivativeQuantityJacobian<ODESystem, Variables, ResidualType>(*odeSystem_, vars, jacobian_);
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

            const auto coordinate = stageParams_->timeAtStage(0);
            Detail::ODE::updateIndependentVariable(
                vars,
                IndependentVariableLevel<Scalar>{coordinate, coordinate, stageParams_->timeStepFraction(0)}
            );

            derivativeQuantities_.push_back(zeroResidualLike_(vars));
            rhsTerms_.push_back(zeroResidualLike_(vars));
            evaluateTerms_(vars, derivativeQuantities_.back(), rhsTerms_.back());
        }

        const auto coordinate = stageParams_->timeAtStage(curStage);
        const auto previousCoordinate = stageParams_->timeAtStage(0);
        const auto stepFraction = stageParams_->timeStepFraction(curStage);
        Detail::ODE::updateIndependentVariable(
            vars,
            IndependentVariableLevel<Scalar>{coordinate, previousCoordinate, stepFraction}
        );

        derivativeQuantities_.push_back(zeroResidualLike_(vars));
        rhsTerms_.push_back(zeroResidualLike_(vars));
    }

    void clearStages()
    {
        derivativeQuantities_.clear();
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
                        ResidualType& derivativeQuantity,
                        ResidualType& rhs) const
    {
        Detail::ODE::setZero(derivativeQuantity);
        Detail::ODE::setZero(rhs);
        Detail::ODE::evaluateDerivativeQuantity(*odeSystem_, vars, derivativeQuantity);
        Detail::ODE::evaluateRhs(*odeSystem_, vars, rhs);
    }

    std::shared_ptr<const ODESystem> odeSystem_;
    ResidualType residual_;
    JacobianMatrix jacobian_;
    SolutionVector previousSolution_;
    std::vector<ResidualType> derivativeQuantities_;
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
    using IntegrationDriver = MultiStageTimeStepper<NonlinearSolver, Scalar>;

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
    , integrationDriver_(nonlinearSolver_, method_, paramGroup)
    {}

    void step(Variables& vars, const Scalar independentVariable, const Scalar stepSize)
    { integrationDriver_.step(vars, independentVariable, stepSize); }

    void solve(Variables& vars, const Scalar endCoordinate, const Scalar stepSize)
    {
        if (!(stepSize > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "ODE integration step size has to be positive");

        auto coordinate = vars.independentVariableLevel().current();
        while (coordinate < endCoordinate)
        {
            using std::min;
            const auto currentStepSize = min(stepSize, endCoordinate - coordinate);
            step(vars, coordinate, currentStepSize);
            coordinate = vars.independentVariableLevel().current();
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
    IntegrationDriver integrationDriver_;
};

} // end namespace Dumux::Experimental

#endif
