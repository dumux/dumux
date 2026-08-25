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
#ifndef DUMUX_ODE_SOLVER_HH
#define DUMUX_ODE_SOLVER_HH

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
 */
template<class Scalar>
class IndependentVariableLevel
{
public:
    /*!
     * \brief Constructs a level at a given value of the independent variable
     *
     * \param current the current value of the independent variable
     */
    explicit IndependentVariableLevel(const Scalar current)
    : current_(current)
    , previous_(current)
    , stepFraction_(1.0)
    {}

    /*!
     * \brief Constructs a level occurring within an integration step
     *
     * \param current      the current value of the independent variable
     * \param previous     the value at the beginning of the integration step
     * \param stepFraction the fraction of the integration step represented by this level
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

/*!
 * \brief Assigns an ODE quantity to an output object of a potentially different type
 *
 * \param  to   destination object
 * \param  from source object
 */
template<class To, class From>
void assign(To& to, const From& from)
{
    if constexpr (std::is_assignable_v<To&, From>)
        to = from;
    else
        static_assert(AlwaysFalse<To>::value, "Cannot assign ODE quantity to the requested output type.");
}

/*!
 * \brief Sets a scalar, vector, or matrix-like value to zero
 *
 * \param value value to set to zero
 */
template<class Value>
void setZero(Value& value)
{
    value = 0.0;
}

/*!
 * \brief Sets a scalar or square matrix-like object to the identity
 *
 * \param matrix object to set to the identity
 */
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

/*!
 * \brief Accumulates the scaled vector \f$y \mathrel{+}= \text{factor} x\f$
 *
 * \param factor scaling factor
 * \param x vector to scale and add
 * \param y vector to update
 */
template<class Vector, class Scalar>
void addScaled(const Scalar factor, const Vector& x, Vector& y)
{
    Dumux::DofBackend<Vector>::axpy(factor, x, y);
}

/*!
 * \brief Scales a matrix-like object in place
 *
 * \param matrix matrix to scale
 * \param factor scaling factor
 */
template<class Matrix, class Scalar>
void scaleMatrix(Matrix& matrix, const Scalar factor)
{
    matrix *= factor;
}

/*!
 * \brief Accumulates the scaled matrix \f$y \mathrel{+}= \text{factor} x\f$
 *
 * \param factor scaling factor
 * \param x      matrix to scale and add
 * \param y      matrix to update
 */
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

/*!
 * \brief Evaluate the right-hand side of an ODE system.
 *
 * \param odeSystem ODE system to evaluate
 * \param vars      variables at which to evaluate the right-hand side
 * \param rhs       right-hand side result
 * \note Supports both `rhs(vars, rhs)` and a value-returning `rhs(vars)`.
 */
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

/*!
 * \brief Evaluate the accumulation mapping \f$M(u,\xi)\f$.
 *
 * \param odeSystem ODE system to evaluate
 * \param vars      variables at which to evaluate the accumulation mapping
 * \param value     resulting accumulation value \f$M(u,\xi)\f$
 * \note If the system provides no `accumulation` function, the primary variables are used,
 *       corresponding to the default mapping \f$M(u,\xi)=u\f$.
 */
template<class ODESystem, class Vars, class Residual>
void evaluateAccumulation(const ODESystem& odeSystem, const Vars& vars, Residual& value)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Residual& result) { system.accumulation(v, result); })
        odeSystem.accumulation(vars, value);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.accumulation(v); })
        assign(value, odeSystem.accumulation(vars));
    else
        assign(value, Dumux::VariablesBackend<Vars>::dofs(vars));
}

/*!
 * \brief Evaluate the Jacobian of the ODE right-hand side.
 *
 * \param odeSystem ODE system to evaluate
 * \param vars      variables at which to evaluate the Jacobian
 * \param jacobian  resulting right-hand-side Jacobian
 * \note Accepts output-argument and value-returning forms named `rhsJacobian`.
 */
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

/*!
 * \brief Evaluate the accumulation Jacobian \f$\partial M/\partial u\f$.
 *
 * \param odeSystem ODE system to evaluate
 * \param vars      variables at which to evaluate the Jacobian
 * \param jacobian  resulting accumulation Jacobian \f$\partial M/\partial u\f$
 * \note Uses the identity for the default mapping \f$M(u,\xi)=u\f$ and requires an
 *       explicit accumulation Jacobian when the ODE system defines a custom accumulation mapping.
 */
template<class ODESystem, class Vars, class Residual, class Jacobian>
void evaluateAccumulationJacobian(const ODESystem& odeSystem, const Vars& vars, Jacobian& jacobian)
{
    if constexpr (requires(const ODESystem& system, const Vars& v, Jacobian& j) { system.accumulationJacobian(v, j); })
        odeSystem.accumulationJacobian(vars, jacobian);
    else if constexpr (requires(const ODESystem& system, const Vars& v) { system.accumulationJacobian(v); })
        assign(jacobian, odeSystem.accumulationJacobian(vars));
    else if constexpr (requires(const ODESystem& system, const Vars& v, Residual& result) { system.accumulation(v, result); }
                       || requires(const ODESystem& system, const Vars& v) { system.accumulation(v); })
        DUNE_THROW(Dune::NotImplemented,
            "ODE systems with a custom accumulation require accumulationJacobian(vars, jacobian) or accumulationJacobian(vars).");
    else
        setIdentity(jacobian);
}

/*!
 * \brief Updates the independent-variable information stored by an ODE variables object
 *
 * \param vars  variables to update
 * \param level new independent-variable level
 */
template<class Vars, class Level>
void updateIndependentVariable(Vars& vars, const Level& level)
{
    if constexpr (requires(Vars& v, const Level& l) { v.updateIndependentVariable(l); })
        vars.updateIndependentVariable(level);
    else
        static_assert(AlwaysFalse<Vars>::value,
            "ODE variables must provide updateIndependentVariable(level) or use Dumux::Experimental::ODEVariables.");
}

/*!
 * \brief Compute the Euclidean norm of a scalar or vector-like object.
 *
 * \param vector the object whose norm is computed
 * \note Uses `two_norm()` when available and otherwise recursively computes the norm.
 */
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

/*!
 * \brief Solve the local linear system \f$A x = b\f$.
 *
 * \param matrix system matrix \f$A\f$
 * \param x      solution vector
 * \param b      right-hand-side vector
 * \return `true` if the supported direct solution operation was performed.
 */
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
    /*!
     * \brief Accepts the Newton residual-reduction setting
     *
     * The residual-reduction argument required by the Newton interface is ignored
     * because the local linear solver needs no adaptation.
     */
    void setResidualReduction(double) {}

    /*!
     * \brief Solve the linearized ODE stage system.
     *
     * \param matrix linearized system matrix
     * \param x      solution vector
     * \param b      right-hand-side vector
     * \return `true` if the local system was solved.
     */
    bool solve(const JacobianMatrix& matrix,
               ResidualVector& x,
               const ResidualVector& b) const
    { return Detail::ODE::solve(matrix, x, b); }

    /*!
     * \brief Computes the norm used by the Newton convergence check
     *
     * \param residual the residual vector
     * \return norm of the residual vector
     */
    auto norm(const ResidualVector& residual) const
    { return Detail::ODE::norm(residual); }
};

/*!
 * \ingroup Experimental
 * \brief Assembler that maps an ODE right-hand-side interface to the DuMux Newton interface.
 *
  * ODE systems are expected in the form
 * \f[
 *     \frac{\mathrm d}{\mathrm d\xi} M(u,\xi) = f(u,\xi),
 * \f]
 * where \f$\xi\f$ denotes the independent variable. If no
 * `accumulation()` is provided, the default is \f$M(u,\xi)=u\f$. The ODE system must export
 * `Scalar`, `SolutionVector`, `ResidualType`, and `JacobianMatrix`. They provide
 * `rhs(vars, rhs)` or `rhs(vars)`. Implicit stages additionally require
 * `rhsJacobian(vars, jacobian)` or `rhsJacobian(vars)`. Optionally, systems can provide
 * `accumulation(vars, value)` and
 * `accumulationJacobian(vars, jacobian)` for equations
 * \f$\frac{\mathrm{d}}{\mathrm{d}\xi} M(u,\xi) = f(u,\xi)\f$, where
 * \f$\xi\f$ denotes an arbitrary independent variable.
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

    /*!
     * \brief Constructs an assembler that shares ownership of an ODE system
     *
     * \param odeSystem the ODE system to assemble
     */
    explicit MultiStageODEAssembler(std::shared_ptr<const ODESystem> odeSystem)
    : odeSystem_(std::move(odeSystem))
    {}

    /*!
     * \brief Constructs an assembler with an internally owned copy of an ODE system
     *
     * \param odeSystem the ODE system to copy
     */
    explicit MultiStageODEAssembler(const ODESystem& odeSystem)
    : MultiStageODEAssembler(std::make_shared<ODESystem>(odeSystem))
    {}

    /*!
     * \brief Satisfies the Newton assembler interface
     *
     * Local ODE systems do not require additional linear-system setup.
     */
    void setLinearSystem() {}

    /*!
     * \brief Returns the assembled Jacobian matrix
     *
     * \return assembled Jacobian matrix
     */
    JacobianMatrix& jacobian()
    { return jacobian_; }

    /*!
     * \brief Returns the assembled Jacobian matrix
     *
     * \return assembled Jacobian matrix
     */
    const JacobianMatrix& jacobian() const
    { return jacobian_; }

    /*!
     * \brief Returns the assembled residual vector
     *
     * \return assembled residual vector
     */
    ResidualType& residual()
    { return residual_; }

    /*!
     * \brief Returns the assembled residual vector
     *
     * \return assembled residual vector
     */
    const ResidualType& residual() const
    { return residual_; }

    /*!
     * \brief Assemble the weighted residual for the currently prepared stage.
     *
     * \param vars variables at which to evaluate the current-stage terms
     */
    void assembleResidual(const Variables& vars)
    {
        if (!stageParams_)
            DUNE_THROW(Dune::InvalidStateException, "No ODE stage prepared before residual assembly");

        if (stageParams_->size() != rhsTerms_.size() || stageParams_->size() != accumulationValues_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of ODE stage residuals");

        evaluateTerms_(vars, accumulationValues_.back(), rhsTerms_.back());

        residual_ = zeroResidualLike_(vars);
        for (std::size_t k = 0; k < stageParams_->size(); ++k)
        {
            if (!stageParams_->skipTemporal(k))
                Detail::ODE::addScaled(stageParams_->temporalWeight(k), accumulationValues_[k], residual_);
            if (!stageParams_->skipSpatial(k))
                Detail::ODE::addScaled(-stageParams_->spatialWeight(k), rhsTerms_[k], residual_);
        }
    }

    /*!
     * \brief Assemble the residual and its Jacobian for the current stage.
     *
     * \param vars variables at which to evaluate the current-stage terms
     */
    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);

        const auto curStage = stageParams_->size() - 1;
        Detail::ODE::evaluateAccumulationJacobian<ODESystem, Variables, ResidualType>(*odeSystem_, vars, jacobian_);
        Detail::ODE::scaleMatrix(jacobian_, stageParams_->temporalWeight(curStage));

        if (!stageParams_->skipSpatial(curStage))
        {
            auto rhsJacobian = jacobian_;
            Detail::ODE::setZero(rhsJacobian);
            Detail::ODE::evaluateRhsJacobian(*odeSystem_, vars, rhsJacobian);
            Detail::ODE::addScaledMatrix(-stageParams_->spatialWeight(curStage), rhsJacobian, jacobian_);
        }
    }

    /*!
     * \brief Prepare the next consecutive stage and update its independent-variable information.
     *
     * \param vars   variables whose independent-variable level is advanced to the stage coordinate
     * \param params coefficients and coordinate data for the stage to prepare
     */
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

            const auto coordinate = stageParams_->timeAtStage(0); // uses existing backend for time integration
            Detail::ODE::updateIndependentVariable(
                vars,
                IndependentVariableLevel<Scalar>{coordinate, coordinate, stageParams_->timeStepFraction(0)}
            );

            accumulationValues_.push_back(zeroResidualLike_(vars));
            rhsTerms_.push_back(zeroResidualLike_(vars));
            evaluateTerms_(vars, accumulationValues_.back(), rhsTerms_.back());
        }

        const auto coordinate = stageParams_->timeAtStage(curStage); // uses existing backend for time integration
        const auto previousCoordinate = stageParams_->timeAtStage(0); // uses existing backend for time integration
        const auto stepFraction = stageParams_->timeStepFraction(curStage);
        Detail::ODE::updateIndependentVariable(
            vars,
            IndependentVariableLevel<Scalar>{coordinate, previousCoordinate, stepFraction}
        );

        accumulationValues_.push_back(zeroResidualLike_(vars));
        rhsTerms_.push_back(zeroResidualLike_(vars));
    }

    /*!
     * \brief Discards all cached stage terms and the active stage parameters
     */
    void clearStages()
    {
        accumulationValues_.clear();
        rhsTerms_.clear();
        stageParams_.reset();
    }

    /*!
     * \brief Returns the solution saved at the beginning of the current integration step
     *
     * \return solution at the beginning of the current integration step
     */
    const SolutionVector& prevSol() const
    { return previousSolution_; }

    /*!
     * \brief Resets stage-local state before repeating or starting an integration step
     *
     * The solution argument required by the integration-driver interface is not needed.
     */
    void resetTimeStep(const SolutionVector&)
    { clearStages(); }

    /*!
     * \brief Returns the ODE system assembled by this object
     *
     * \return assembled ODE system
     */
    const ODESystem& odeSystem() const
    { return *odeSystem_; }

private:
    using Backend = Dumux::VariablesBackend<Variables>;
    using ResidualBackend = Dumux::DofBackend<ResidualType>;

    /*!
     * \brief Creates a zero residual with a size compatible with the given variables
     *
     * \param vars the variables defining the required residual size
     * \return a zero-initialized residual
     */
    ResidualType zeroResidualLike_(const Variables& vars) const
    {
        auto result = ResidualBackend::zeros(Backend::size(Backend::dofs(vars)));
        Detail::ODE::setZero(result);
        return result;
    }

    /*!
     * \brief Evaluates the accumulation and right-hand-side values at the given variables
     *
     * \param vars         variables at which to evaluate the values
     * \param accumulation resulting accumulation value \f$M(u,\xi)\f$
     * \param rhs          resulting right-hand-side value \f$f(u,\xi)\f$
     */
    void evaluateTerms_(const Variables& vars,
                        ResidualType& accumulation,
                        ResidualType& rhs) const
    {
        Detail::ODE::setZero(accumulation);
        Detail::ODE::setZero(rhs);
        Detail::ODE::evaluateAccumulation(*odeSystem_, vars, accumulation);
        Detail::ODE::evaluateRhs(*odeSystem_, vars, rhs);
    }

    std::shared_ptr<const ODESystem> odeSystem_;
    ResidualType residual_;
    JacobianMatrix jacobian_;
    SolutionVector previousSolution_;
    std::vector<ResidualType> accumulationValues_;
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

    /*!
     * \brief Construct a solver with a default linear solver and MPI communication.
     *
     * \param odeSystem  ODE system to solve.
     * \param method     multi-stage integration method.
     * \param paramGroup parameter group used to configure the solver components
     */
    ODESolver(std::shared_ptr<const ODESystem> odeSystem,
              std::shared_ptr<const Method> method,
              const std::string& paramGroup = "")
    : ODESolver(std::move(odeSystem),
                std::make_shared<LinearSolver>(),
                std::move(method),
                Dune::MPIHelper::getCommunication(),
                paramGroup)
    {}

    /*!
     * \brief Construct a solver with a default linear solver and custom communication.
     *
     * \param odeSystem  ODE system to solve
     * \param method     multi-stage integration method
     * \param comm       communication object used by the nonlinear solver
     * \param paramGroup parameter group used to configure the solver components
     */
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

    /*!
     * \brief Construct a solver with a custom linear solver and MPI communication.
     *
     * \param odeSystem    ODE system to solve
     * \param linearSolver linear solver used for implicit stages
     * \param method       multi-stage integration method
     * \param paramGroup   parameter group used to configure the solver components
     */
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

    /*!
     * \brief Construct a solver with custom linear solver and communication objects.
     *
     * \param odeSystem    ODE system to solve
     * \param linearSolver linear solver used for implicit stages
     * \param method       multi-stage integration method
     * \param comm         communication object used by the nonlinear solver
     * \param paramGroup   parameter group used to configure the solver components
     */
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

    /*!
     * \brief Advance the variables by one integration step.
     *
     * \param vars                variables to update in place
     * \param independentVariable value of the independent variable at the beginning of the step
     * \param stepSize            integration-step size
     */
    void step(Variables& vars, const Scalar independentVariable, const Scalar stepSize)
    { integrationDriver_.step(vars, independentVariable, stepSize); }

    /*!
     * \brief Integrate repeatedly with a fixed maximum step size up to an end coordinate.
     *
     * \param vars          variables to update in place
     * \param endCoordinate final value of the independent variable
     * \param stepSize      requested step size; the last step may be shortened
     */
    void solve(Variables& vars, const Scalar endCoordinate, const Scalar stepSize)
    {
        if (!(stepSize > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "ODE integration step size has to be positive");

        auto startCoordinate = vars.independentVariableLevel().current();

        if (endCoordinate < startCoordinate)
            DUNE_THROW(
                Dune::InvalidStateException,
                "Backward ODE integration is not supported: end coordinate "
                << endCoordinate << " is smaller than start coordinate "
                << startCoordinate
            );

        if (endCoordinate == startCoordinate)
            return;

        while (startCoordinate < endCoordinate)
        {
            using std::min;
            const auto currentStepSize = min(stepSize, endCoordinate - startCoordinate);
            step(vars, startCoordinate, currentStepSize);
            startCoordinate = vars.independentVariableLevel().current();
        }
    }

    /*!
     * \brief Returns the ODE assembler
     *
     * \return ODE assembler
     */
    Assembler& assembler()
    { return *assembler_; }

    /*!
     * \brief Returns the ODE assembler
     *
     * \return ODE assembler
     */
    const Assembler& assembler() const
    { return *assembler_; }

    /*!
     * \brief Returns the nonlinear solver used for implicit stages
     *
     * \return nonlinear solver
     */
    NonlinearSolver& nonlinearSolver()
    { return *nonlinearSolver_; }

    /*!
     * \brief Returns the nonlinear solver used for implicit stages
     *
     * \return nonlinear solver
     */
    const NonlinearSolver& nonlinearSolver() const
    { return *nonlinearSolver_; }

    /*!
     * \brief Returns the linear solver used by the nonlinear solver
     *
     * \return linear solver
     */
    LinearSolver& linearSolver()
    { return *linearSolver_; }

    /*!
     * \brief Returns the linear solver used by the nonlinear solver
     *
     * \return linear solver
     */
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
