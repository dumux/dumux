// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Newmark-beta time integration scheme
 */

#ifndef DUMUX_TIMESTEPPING_NEWMARK_BETA_HH
#define DUMUX_TIMESTEPPING_NEWMARK_BETA_HH

#include <dumux/common/parameters.hh>

namespace Dumux::Experimental {

/*!
 * \brief Newmark-beta time integration scheme
 * \tparam Scalar The scalar type used for the solution
 * \tparam Solution The solution type used for the variables
 * \note Newmark (1959) "A method of computation for structural dynamics"
 *       Journal of the Engineering Mechanics Division,
 *       https://doi.org/10.1061/JMCEA3.0000098
 * \note This is typically used for PDEs of the form M*a + C*v + K*x = f
 *       with M the mass matrix, C the damping matrix, K the stiffness matrix,
 *       x the displacement, v the velocity and a the acceleration.
 *       An example is dynamic mechanics.
 */
template<class Scalar, class Solution>
class NewmarkBeta
{
    using Block = typename Solution::block_type;
public:
    //! Construct the Newmark-beta time integration scheme
    //! with the given β and γ parameters
    NewmarkBeta(const Scalar beta, const Scalar gamma)
    : beta_(beta)
    , gamma_(gamma)
    {}

    //! Construct the Newmark-beta time integration scheme
    //! where parameters are read from the parameter file
    NewmarkBeta()
    : NewmarkBeta(
        getParam<Scalar>("NewmarkBeta.Beta", 0.25),
        getParam<Scalar>("NewmarkBeta.Gamma", 0.5)
    ) {}

    //! Initialize the time integration scheme with the current solution
    //! while setting the velocity and acceleration to zero
    void initialize(const Solution& x)
    {
        x_ = x;

        // make sure the size is correct, but the values are zero
        v_ = x; v_ = 0.0;
        a_ = x; a_ = 0.0;
    }

    //! Initialize the time integration scheme with the current solution
    void initialize(const Solution& x, const Solution& v, const Solution& a)
    {
        x_ = x;
        v_ = v;
        a_ = a;
    }

    //! Update with new position x
    void update(const Scalar dt, const Solution& x)
    {
        for (std::size_t i = 0; i < x_.size(); ++i)
        {
            const auto xNew = x[i];
            const auto aNew = acceleration(i, dt, xNew);
            const auto vNew = velocity(i, dt, xNew, aNew);

            x_[i] = xNew;
            v_[i] = vNew;
            a_[i] = aNew;
        }
    }

    //! new a in terms of the old x, v, a, the new x and the time step size
    Block acceleration(const std::size_t index, const Scalar dt, const Block& x) const
    {
        const auto& x0 = x_[index];
        const auto& v0 = v_[index];
        const auto& a0 = a_[index];

        return (x - x0 - dt*v0)/(beta_*dt*dt) + (beta_ - 0.5)/beta_*a0;
    }

    //! new v in terms of the old v, a, the new x, a and the time step size
    Block velocity(const std::size_t index, const Scalar dt, const Block& x, const Block& a) const
    {
        const auto& v0 = v_[index];
        const auto& a0 = a_[index];

        return v0 + dt*((1.0-gamma_)*a0 + gamma_*a);
    }

    //! new v in terms of the old v, a, the new x and the time step size
    Block velocity(const std::size_t index, const Scalar dt, const Block& x) const
    {
        const auto a = acceleration(index, dt, x);
        return velocity(index, dt, x, a);
    }

    //! current position
    Scalar position(const std::size_t index) const
    { return x_[index]; }

    //! current velocity
    Scalar velocity(const std::size_t index) const
    { return v_[index]; }

    //! current acceleration
    Scalar acceleration(const std::size_t index) const
    { return a_[index]; }

private:
    Scalar beta_, gamma_; //!< Newmark-beta parameters
    Solution x_, v_, a_; //!< current position, velocity and acceleration
};

} // end namespace Dumux::Experimental

#endif
