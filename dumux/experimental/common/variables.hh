// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Core
 * \copydoc Dumux::Experimental::Variables
 */
#ifndef DUMUX_COMMON_VARIABLES_HH
#define DUMUX_COMMON_VARIABLES_HH

#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dumux/experimental/timestepping/timelevel.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Core
 * \brief Class that represents the variables of a model.
 *        We assume that models are formulated on the basis of primary and
 *        possibly secondary variables, where the latter may non-linearly
 *        depend on the former. Variables in Dumux represent the state of
 *        a numerical solution of a model, consisting of all primary/secondary
 *        variables and, if a transient problem is modeled, time information.
 *
 *        This class defines the interface that is expected of variable classes,
 *        and it provides the implementation for models that do not require storing
 *        any additional information besides the primary variables and (optionally)
 *        time.
 * \tparam X The type used for solution vectors, i.e. all primary variables.
 */
template<class X>
class Variables
{
    template<class T, bool indexable = Dune::IsIndexable<T>::value>
    struct ScalarT { using type = T; };

    template<class T>
    struct ScalarT<T, true>
    {
        using Element = std::decay_t<decltype(std::declval<T>()[0])>;
        using type = typename ScalarT<Element>::type;
    };

public:
    //! export the type of solution vector
    using SolutionVector = X;

    //! export the underlying scalar type
    using Scalar = typename ScalarT<X>::type;

    //! export the time representation
    using TimeLevel = Dumux::Experimental::TimeLevel<Scalar>;

    //! Default constructor
    explicit Variables() : x_(), t_(0.0) {}

    //! Construction from a solution
    explicit Variables(const SolutionVector& x,
                       const TimeLevel& t = TimeLevel{0.0})
    : x_(x), t_(t)
    {}

    //! Construction from a solution
    explicit Variables(SolutionVector&& x,
                       const TimeLevel& t = TimeLevel{0.0})
    : x_(std::move(x)), t_(t)
    {}

    //! Construction from initializer lambda
    template<class Initializer,
              std::enable_if_t<(std::is_invocable_r_v<void, Initializer, X&>), int> = 0>
    explicit Variables(const Initializer& initializeSolution,
                       const TimeLevel& timeLevel = TimeLevel{0.0})
    : t_(timeLevel)
    {
        initializeSolution(x_);
    }

    //! Return the time level
    const TimeLevel& timeLevel() const
    { return t_; }

    //! Return reference to the solution
    const SolutionVector& dofs() const { return x_; }

    //! Non-const access still required for privar switch (TODO: Remove dependency)
    SolutionVector& dofs() { return x_; }

    //! Update the state to a new solution
    void update(const SolutionVector& x)
    { x_ = x; }

    //! Update the time level only
    void updateTime(const TimeLevel& t)
    { t_ = t; }

    //! Update the state to a new solution & time level
    void update(const SolutionVector& x,
                const TimeLevel& t)
    {
        x_ = x;
        t_ = t;
    }

private:
    SolutionVector x_;
    TimeLevel t_;
};

} // end namespace Dumux::Experimental

#endif
