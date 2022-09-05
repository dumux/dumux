// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Common
 * \brief Default implementation for basic variables.
 */
#ifndef DUMUX_COMMON_DEFAULT_VARIABLES_HH
#define DUMUX_COMMON_DEFAULT_VARIABLES_HH

#include <dumux/experimental/new_assembly/dumux/common/timelevel.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Default implementation for basic (possibly time-dependent) variables.
 * \tparam X The type used for storing the degrees of freedom. Must be a type
 *           that is interoperable with our linear system functions.
 */
template<LinearSystem::Interoperable X>
class DefaultVariables
{
    using Scalar = LinearSystem::ScalarType<X>;

public:
    using Dofs = X;
    using TimeLevel = Dumux::TimeLevel<Scalar>;

    explicit DefaultVariables()
    : x_(), t_(0.0)
    {}

    explicit DefaultVariables(const Dofs& x,
                              const TimeLevel& t = TimeLevel{0.0})
    : x_(x), t_(t)
    {}

    explicit DefaultVariables(Dofs&& x,
                              const TimeLevel& t = TimeLevel{0.0})
    : x_(std::move(x)), t_(t)
    {}

    template<std::invocable<Dofs&> Initializer>
    explicit DefaultVariables(const Initializer& initializeSolution,
                              const TimeLevel& timeLevel = TimeLevel{0.0})
    : t_(timeLevel)
    { initializeSolution(x_); }

    //! Return the time level
    const TimeLevel& timeLevel() const
    { return t_; }

    //! Return reference to the solution
    const Dofs& dofs() const
    { return x_; }

    //! Update the state to a new solution
    void update(const Dofs& x)
    { x_ = x; }

    //! Update the time level only
    void updateTime(const TimeLevel& t)
    { t_ = t; }

    //! Update the state to a new solution & time level
    void update(const Dofs& x, const TimeLevel& t)
    {
        x_ = x;
        t_ = t;
    }

private:
    Dofs x_;
    TimeLevel t_;
};

} // namespace Dumux

#endif
