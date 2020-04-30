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
 * \brief Parameters for different multistage time stepping methods
 * \note See e.g. https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_METHODS_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_METHODS_HH

#include <cmath>
#include <string>
#include <array>

namespace Dumux {

/*!
 * \brief Abstract interface for one-step multi-stage method parameters in Shu/Osher form.
 *
 * This implementation is based on the Shu/Osher form from:
 * Chi W. Shu and Stanley Osher. Efficient implementation of essentially
 * non- oscillatory shock-capturing schemes. J. Comput. Phys., 77:439–471, 1988.
 * https://doi.org/10.1016/0021-9991(88)90177-5.
 * To this end Eq. (2.12) is extended for implicit schemes.
 *
 * We consider the general PDE form
 *
 * \f[
 * \begin{equation}
 *   \frac{\partial M(x)}{\partial t} - R(x, t) = 0,
 * \end{equation}
 * \f]
 *
 * where \f$ M(x)\f$ is the temporal operator/residual in terms of the solution \f$ x \f$,
 * and \f$ R(x)\f$ is the discrete spatial operator/residual.
 * \f$ M(x)\f$ usually corresponds to the conserved quanity (e.g. mass), and \f$ R(x)\f$
 * contains the rest of the residual. We can then construct \f$ m \f$-stage time discretization methods.
 * Integrating from time \f$ t^n\f$ to \f$ t^{n+1}\f$ with time step size \f$ \Delta t^n\f$, we solve:
 *
 * \f[
 * \begin{aligned}
 *   x^{(0)} &= u^n\\
 *   \sum_{k=0}^i \left[ \alpha_{ik} M\left(x^{(k)}, t^n + d_k\Delta t^n\right)
 *     + \beta_{ik}\Delta t^n R \left(x^{(k)}, t^n+d_k\Delta t^n \right)\right] &= 0 & \forall i \in \{1,\ldots,m\} \\
 *   x^{n+1} &= x^{(m)}
 * \end{aligned}
 * \f]
 * where \f$ x^{(k)} \f$ denotes the intermediate solution at stage \f$ k \f$.
 * Dependent on the number of stages \f$ m \f$, and the coefficients \f$ \alpha, \beta, d\f$,
 * schemes with different properties and order of accuracy can be constructed.
 */
template<class Scalar>
class MultiStageMethod
{
public:
    virtual bool implicit () const = 0;

    virtual std::size_t numStages () const = 0;

    //! weights of the temporal operator residual (\f$ \alpha_{ik} \f$)
    virtual Scalar temporalWeight (std::size_t i, std::size_t k) const = 0;

    //! weights of the spatial operator residual (\f$ \beta_{ik} \f$)
    virtual Scalar spatialWeight (std::size_t i, std::size_t k) const = 0;

    //! time step weights for each stage (\f$ d_k \f$)
    virtual Scalar timeStepWeight (std::size_t k) const = 0;

    virtual std::string name () const = 0;

    virtual ~MultiStageMethod() = default;
};

//! Multi-stage time stepping scheme implementations
namespace MultiStage {

/*!
 * \brief A theta time stepping scheme
 * theta=1.0 is an implicit Euler scheme,
 * theta=0.0 an explicit Euler scheme,
 * theta=0.5 is a Cranck-Nicholson scheme
 */
template<class Scalar>
class Theta : public MultiStageMethod<Scalar>
{
public:
    explicit Theta(const Scalar theta)
    : paramAlpha_{{-1.0, 1.0}}
    , paramBeta_{{1.0-theta, theta}}
    , paramD_{{0.0, 1.0}};
    {}

    bool implicit () const final
    { return paramBeta_[1] > 0.0; }

    std::size_t numStages () const final
    { return 1; }

    Scalar temporalWeight (std::size_t, std::size_t k) const final
    { return paramAlpha_[k]; }

    Scalar spatialWeight (std::size_t, std::size_t k) const final
    { return paramBeta_[k]; }

    Scalar timeStepWeight (std::size_t k) const final
    { return paramD_[k]; }

    std::string name () const override
    { return "theta scheme"; }

private:
    std::array<Scalar, 2> paramAlpha_;
    std::array<Scalar, 2> paramBeta_;
    std::array<Scalar, 2> paramD_;
};

/*!
 * \brief An explicit Euler time stepping scheme
 */
template<class Scalar>
class ExplicitEuler final : public Theta<Scalar>
{
public:
    ExplicitEuler() : Theta<Scalar>(0.0) {}

    std::string name () const final
    { return "explicit Euler"; }
};

/*!
 * \brief An implicit Euler time stepping scheme
 */
template<class Scalar>
class ImplicitEuler final : public Theta<Scalar>
{
public:
    ImplicitEuler() : Theta<Scalar>(1.0) {}

    std::string name () const final
    { return "implicit Euler"; }
};

/*!
 * \brief Classical explicit fourth order Runge-Kutta scheme
 */
template<class Scalar>
class RungeKuttaExplicitFourthOrder final : public MultiStageMethod<Scalar>
{
    RungeKuttaExplicitFourthOrder()
    : paramAlpha_{{{-1.0, 1.0, 0.0, 0.0, 0.0},
                   {-1.0, 0.0, 1.0, 0.0, 0.0},
                   {-1.0, 0.0, 0.0, 1.0, 0.0},
                   {-1.0, 0.0, 0.0, 0.0, 1.0}}}
    , paramBeta_{{{0.5, 0.0, 0.0, 0.0, 0.0},
                  {0.0, 0.5, 0.0, 0.0, 0.0},
                  {0.0, 0.0, 1.0, 0.0, 0.0},
                  {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0.0}}}
    , paramD_{{0.0, 0.5, 0.5, 1.0, 1.0}};
    {}

    bool implicit () const final
    { return false; }

    std::size_t numStages () const final
    { return 4; }

    Scalar temporalWeight (std::size_t i, std::size_t k) const final
    { return paramAlpha_[i-1][k]; }

    Scalar spatialWeight (std::size_t i, std::size_t k) const final
    { return paramBeta_[i-1][k]; }

    Scalar timeStepWeight (std::size_t k) const final
    { return paramD_[k]; }

    std::string name () const final
    { return "explicit Runge-Kutta 4th order"; }

private:
    std::array<std::array<Scalar, 5>, 4> paramAlpha_;
    std::array<std::array<Scalar, 5>, 4> paramBeta_;
    std::array<Scalar, 5> paramD_;
};

} // end namespace MultiStage
} // end namespace Dumux

#endif
