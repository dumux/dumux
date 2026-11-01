// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Data exchanged between time steppers and assemblers
 *
 * We consider the semi-discrete problem
 *
 * \f[
 *   \frac{\partial^2 I(x, t)}{\partial t^2} + \frac{\partial M(x, t)}{\partial t} + R(x, t) = 0,
 * \f]
 *
 * with the inertia operator \f$ I \f$ (zero unless the problem enables it), the storage
 * operator \f$ M \f$ and the spatial operator \f$ R \f$ (fluxes and sources). A time stepper
 * reduces each time step to a sequence of implicit solves of the form
 *
 * \f[
 *   w_I I(x, t) + w_M M(x, t) + w_R R(x, t) + c = 0,
 * \f]
 *
 * where the vector \f$ c \f$ collects all contributions of states that are
 * already known (previous stages or previous time steps). The assembler
 * evaluates this residual and its Jacobian
 * \f$ w_I \partial I / \partial x + w_M \partial M / \partial x + w_R \partial R / \partial x \f$.
 * It holds no time integration history: the stepper evaluates the
 * unweighted operators at the states it needs to remember and builds
 * \f$ c \f$ from them.
 *
 * Time-dependent problems receive the time at which the assembler evaluates
 * the operators through `void setTime(Scalar t) const`.
 */
#ifndef DUMUX_TIMESTEPPING_IMPLICIT_SOLVE_HH
#define DUMUX_TIMESTEPPING_IMPLICIT_SOLVE_HH

#include <memory>
#include <concepts>
#include <type_traits>

#include <dumux/experimental/timestepping/timelevel.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \brief Weights of the operators in the residual of an implicit solve
 * \note An operator with a weight of exactly zero is not evaluated.
 * \note The inertia operator \f$ I \f$ enters problems with second time derivatives,
 *       \f$ \partial^2 I(x, t)/\partial t^2 + \partial M(x, t)/\partial t + R(x, t) = 0 \f$.
 *       It is evaluated only for problems that enable it (`static constexpr bool enableInertia()`).
 */
template<class Scalar>
struct OperatorWeights
{
    Scalar inertia = 0.0; //!< weight \f$ w_I \f$ of the inertia operator
    Scalar storage = 0.0; //!< weight \f$ w_M \f$ of the storage operator
    Scalar spatial = 0.0; //!< weight \f$ w_R \f$ of the spatial operator
};

/*!
 * \ingroup Experimental
 * \brief Everything an assembler needs to know about an implicit solve
 */
template<class Scalar, class Residual>
struct ImplicitSolveParams
{
    //! the time level at which the operators are evaluated
    TimeLevel<Scalar> timeLevel;

    //! the weights of the operators
    OperatorWeights<Scalar> weights;

    //! the contribution \f$ c \f$ of known states (nullptr means zero)
    std::shared_ptr<const Residual> constant;
};

/*!
 * \ingroup Experimental
 * \brief Output vectors for the unweighted operator evaluations (nullptr means not requested)
 */
template<class Residual>
struct OperatorEvaluations
{
    Residual* inertia = nullptr;
    Residual* storage = nullptr;
    Residual* spatial = nullptr;
};

} // end namespace Dumux::Experimental

namespace Dumux::Experimental::Detail::TimeStepping {

/*!
 * \brief Whether the problem enables the inertia operator
 * \note A problem enables inertia with `static constexpr bool enableInertia()`.
 */
template<class Problem>
constexpr bool enablesInertia()
{
    if constexpr (requires { { Problem::enableInertia() } -> std::convertible_to<bool>; })
        return Problem::enableInertia();
    else
        return false;
}

struct TimeSetterNames
{
    void setTime();
    void updateTime();
    void setTimeLoop();
};

template<class Problem>
struct ProblemWithTimeSetterNames : public Problem, public TimeSetterNames {};

/*!
 * \brief Whether the problem has a member named setTime, updateTime or setTimeLoop (with any signature)
 * \note Such a member makes the name lookup in ProblemWithTimeSetterNames ambiguous.
 *       Final classes cannot be derived from and are reported as having no such member.
 */
template<class Problem>
constexpr bool hasTimeSetter()
{
    if constexpr (std::is_class_v<Problem> && !std::is_final_v<Problem>)
    {
        using P = ProblemWithTimeSetterNames<Problem>;
        return !requires { &P::setTime; }
            || !requires { &P::updateTime; }
            || !requires { &P::setTimeLoop; };
    }
    else
        return false;
}

/*!
 * \brief Set the time at which the operators of the problem are evaluated
 * \note Time-dependent problems receive the time through `void setTime(Scalar t) const`.
 *       A problem without a member named setTime, updateTime or setTimeLoop is time-independent.
 */
template<class Problem, class Scalar>
void setProblemTime(const Problem& problem, const Scalar t)
{
    if constexpr (requires { problem.setTime(t); })
        problem.setTime(t);
    else
        static_assert(!hasTimeSetter<Problem>(),
            "Time-dependent problems receive the time at which the operators are evaluated through void setTime(Scalar t) const");
}

//! y += a*x for scalars and (block) vectors
template<class Scalar, class Vector>
void axpy(Vector& y, const Scalar a, const Vector& x)
{
    if constexpr (requires { y.axpy(a, x); })
        y.axpy(a, x);
    else
        y += a*x;
}

} // end namespace Dumux::Experimental::Detail::TimeStepping

#endif
