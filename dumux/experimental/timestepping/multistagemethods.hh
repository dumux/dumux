// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright © dune-pdelab developers, see LICENSE.md (permalink below)
// SPDX-License-Identifier: GPL-3.0-or-later
// The code is based on the implementation of time stepping method parameters
// in dune-pdelab (https://archive.softwareheritage.org/swh:1:cnt:9c3d72412f8e4d48d84a0090b2bb4362b5a0d843)
// licensed under GPL-3.0-or-later, see their LICENSE.md for a full list of copyright holders at
// https://archive.softwareheritage.org/swh:1:cnt:b11b484e74eefe20c74f7043309b1b02853df2eb.
// Modifications (different interface naming, comments, types) are licensed under GPL-3.0-or-later.
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Parameters for different multistage time stepping methods
 * \note See e.g. https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_METHODS_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_METHODS_HH

#include <cmath>
#include <string>
#include <array>

namespace Dumux::Experimental {

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
 *   \frac{\partial M(x, t)}{\partial t} - R(x, t) = 0,
 * \end{equation}
 * \f]
 *
 * where \f$ M(x, t)\f$ is the temporal operator/residual in terms of the solution \f$ x \f$,
 * and \f$ R(x, t)\f$ is the discrete spatial operator/residual.
 * \f$ M(x)\f$ usually corresponds to the conserved quantity (e.g. mass), and \f$ R(x)\f$
 * contains the rest of the residual. We can then construct \f$ m \f$-stage time discretization methods.
 * Integrating from time \f$ t^n\f$ to \f$ t^{n+1}\f$ with time step size \f$ \Delta t^n\f$, we solve:
 *
 * \f[
 * \begin{aligned}
 *   x^{(0)} &= u^n\\
 *   \sum_{k=0}^i \left[ \alpha_{ik} M\left(x^{(k)}, t^n + d_k\Delta t^n\right)
 *     + \beta_{ik}\Delta t^n R \left(x^{(k)}, t^n+d_k\Delta t^n \right)\right]
 *      &= 0 & \forall i \in \{1,\ldots,m\} \\
 *   x^{n+1} &= x^{(m)}
 * \end{aligned}
 * \f]
 * where \f$ x^{(k)} \f$ denotes the intermediate solution at stage \f$ k \f$.
 * Dependent on the number of stages \f$ m \f$, and the coefficients \f$ \alpha, \beta, d\f$,
 * schemes with different properties and order of accuracy can be constructed.
 *
 * That the summation only goes up to \f$ i \f$ in stage \f$ i \f$ means that we
 * restrict ourselves to diagonally-implicit Runge-Kutta schemes (DIRK)
 * and explicit schemes.
 *
 * Note that when computing the Jacobian of the residual with respect
 * to \f$ x^{(k)} \f$ at stage \f$ k \f$, only the terms containing the solution of
 * the current stage \f$ k \f$ contribute to the derivatives.
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
    , paramD_{{0.0, 1.0}}
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
 * \brief Second order explicit Runge-Kutta scheme (Heun / explicit trapezoidal rule)
 */
template<class Scalar>
class RungeKuttaExplicitSecondOrderHeun final : public MultiStageMethod<Scalar>
{
public:
    RungeKuttaExplicitSecondOrderHeun()
    : paramAlpha_{{{-1.0, 1.0, 0.0},
                   {-1.0, 0.0, 1.0}}}
    , paramBeta_{{{1.0, 0.0, 0.0},
                  {0.5, 0.5, 0.0}}}
    , paramD_{{0.0, 1.0, 1.0}}
    {}

    bool implicit () const final
    { return false; }

    std::size_t numStages () const final
    { return 2; }

    Scalar temporalWeight (std::size_t i, std::size_t k) const final
    { return paramAlpha_[i-1][k]; }

    Scalar spatialWeight (std::size_t i, std::size_t k) const final
    { return paramBeta_[i-1][k]; }

    Scalar timeStepWeight (std::size_t k) const final
    { return paramD_[k]; }

    std::string name () const final
    { return "explicit Runge-Kutta 2nd order (Heun)"; }

private:
    std::array<std::array<Scalar, 3>, 2> paramAlpha_;
    std::array<std::array<Scalar, 3>, 2> paramBeta_;
    std::array<Scalar, 3> paramD_;
};

/*!
 * \brief Classical explicit third order Runge-Kutta scheme
 */
template<class Scalar>
class RungeKuttaExplicitThirdOrder final : public MultiStageMethod<Scalar>
{
public:
    RungeKuttaExplicitThirdOrder()
    : paramAlpha_{{{-1.0, 1.0, 0.0, 0.0},
                   {-1.0, 0.0, 1.0, 0.0},
                   {-1.0, 0.0, 0.0, 1.0}}}
    , paramBeta_{{{0.5, 0.0, 0.0, 0.0},
                  {-1.0, 2.0, 0.0, 0.0},
                  {1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0}}}
    , paramD_{{0.0, 0.5, 1.0, 1.0}}
    {}

    bool implicit () const final
    { return false; }

    std::size_t numStages () const final
    { return 3; }

    Scalar temporalWeight (std::size_t i, std::size_t k) const final
    { return paramAlpha_[i-1][k]; }

    Scalar spatialWeight (std::size_t i, std::size_t k) const final
    { return paramBeta_[i-1][k]; }

    Scalar timeStepWeight (std::size_t k) const final
    { return paramD_[k]; }

    std::string name () const final
    { return "explicit Runge-Kutta 3rd order"; }

private:
    std::array<std::array<Scalar, 4>, 3> paramAlpha_;
    std::array<std::array<Scalar, 4>, 3> paramBeta_;
    std::array<Scalar, 4> paramD_;
};

/*!
 * \brief Classical explicit fourth order Runge-Kutta scheme
 */
template<class Scalar>
class RungeKuttaExplicitFourthOrder final : public MultiStageMethod<Scalar>
{
public:
    RungeKuttaExplicitFourthOrder()
    : paramAlpha_{{{-1.0, 1.0, 0.0, 0.0, 0.0},
                   {-1.0, 0.0, 1.0, 0.0, 0.0},
                   {-1.0, 0.0, 0.0, 1.0, 0.0},
                   {-1.0, 0.0, 0.0, 0.0, 1.0}}}
    , paramBeta_{{{0.5, 0.0, 0.0, 0.0, 0.0},
                  {0.0, 0.5, 0.0, 0.0, 0.0},
                  {0.0, 0.0, 1.0, 0.0, 0.0},
                  {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0.0}}}
    , paramD_{{0.0, 0.5, 0.5, 1.0, 1.0}}
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


/*!
 * \brief Second order singly diagonally implicit Runge-Kutta scheme
 * \note See Alexander (1977) https://doi.org/10.1137/0714068.
 */
template<class Scalar>
class DIRKSecondOrderAlexander final : public MultiStageMethod<Scalar>
{
public:
    DIRKSecondOrderAlexander()
    {
        using std::sqrt;
        const Scalar gamma = 1.0 - 1.0/sqrt(2.0);

        paramD_ = {{0.0, gamma, 1.0}};
        paramAlpha_ = {{
            {-1.0, 1.0, 0.0},
            {-1.0, 0.0, 1.0}
        }};
        paramBeta_ = {{
            {0.0, gamma, 0.0},
            {0.0, 1.0-gamma, gamma}
        }};
    }

    bool implicit () const final
    { return true; }

    std::size_t numStages () const final
    { return 2; }

    Scalar temporalWeight (std::size_t i, std::size_t k) const final
    { return paramAlpha_[i-1][k]; }

    Scalar spatialWeight (std::size_t i, std::size_t k) const final
    { return paramBeta_[i-1][k]; }

    Scalar timeStepWeight (std::size_t k) const final
    { return paramD_[k]; }

    std::string name () const final
    { return "diagonally implicit Runge-Kutta 2nd order (Alexander)"; }

private:
    std::array<std::array<Scalar, 3>, 2> paramAlpha_;
    std::array<std::array<Scalar, 3>, 2> paramBeta_;
    std::array<Scalar, 3> paramD_;
};

/*!
 * \brief Third order DIRK scheme
 * \note see Alexander (1977) https://doi.org/10.1137/0714068 (Theorem 5)
 */
template<class Scalar>
class DIRKThirdOrderAlexander final : public MultiStageMethod<Scalar>
{
public:
    DIRKThirdOrderAlexander()
    {
        constexpr Scalar alpha = []{
            // Newton iteration for alpha
            Scalar alpha = 0.4358665215; // initial guess
            for (int i = 0; i < 10; ++i)
                alpha = alpha - (alpha*(alpha*alpha-3.0*(alpha-0.5))-1.0/6.0)/(3.0*alpha*(alpha-2.0)+1.5);
            return alpha;
        }();

        constexpr Scalar tau2 = (1.0+alpha)*0.5;
        constexpr Scalar b1 = -(6.0*alpha*alpha -16.0*alpha + 1.0)*0.25;
        constexpr Scalar b2 = (6*alpha*alpha - 20.0*alpha + 5.0)*0.25;

        paramD_ = {{0.0, alpha, tau2, 1.0}};
        paramAlpha_ = {{
            {-1.0, 1.0, 0.0, 0.0},
            {-1.0, 0.0, 1.0, 0.0},
            {-1.0, 0.0, 0.0, 1.0}
        }};
        paramBeta_ = {{
            {0.0, alpha, 0.0, 0.0},
            {0.0, tau2-alpha, alpha, 0.0},
            {0.0, b1, b2, alpha}
        }};
    }

    bool implicit () const final
    { return true; }

    std::size_t numStages () const final
    { return 3; }

    Scalar temporalWeight (std::size_t i, std::size_t k) const final
    { return paramAlpha_[i-1][k]; }

    Scalar spatialWeight (std::size_t i, std::size_t k) const final
    { return paramBeta_[i-1][k]; }

    Scalar timeStepWeight (std::size_t k) const final
    { return paramD_[k]; }

    std::string name () const final
    { return "diagonally implicit Runge-Kutta 3rd order (Alexander)"; }

private:
    std::array<std::array<Scalar, 4>, 3> paramAlpha_;
    std::array<std::array<Scalar, 4>, 3> paramBeta_;
    std::array<Scalar, 4> paramD_;
};

/*!
 * \brief Second order symplectic diagonally implicit Runge-Kutta scheme (Qin-Zhang)
 * \note See Qin & Zhang (1992) \cite QinZhang1992, the two-stage second-order
 *       symplectic DIRK. Being symplectic it satisfies the Cooper condition
 *       \f$ b_i a_{ij} + b_j a_{ji} - b_i b_j = 0 \f$ and therefore conserves all
 *       quadratic invariants exactly (e.g. angular momentum and the energy of a
 *       linear Hamiltonian system). It is A-stable but not L-stable
 *       (\f$ |R(z)| \to 1 \f$ as \f$ \mathrm{Re}(z) \to -\infty \f$), so it does
 *       not damp under-resolved high-frequency modes.
 *
 * Butcher tableau (c | A ; b):
 * \verbatim
 *   1/4 | 1/4   0
 *   3/4 | 1/2  1/4
 *   ----+----------
 *       | 1/2  1/2
 * \endverbatim
 * The method is not stiffly accurate (\f$ a_{sj} \neq b_j \f$), so the final
 * update \f$ x^{n+1} = x^n + \Delta t\,(b_1 R_1 + b_2 R_2) \f$ is appended as an
 * explicit third stage (\f$ \beta_{33} = 0 \f$).
 */
template<class Scalar>
class QinZhangSymplecticDIRK final : public MultiStageMethod<Scalar>
{
public:
    QinZhangSymplecticDIRK()
    : paramAlpha_{{
        {-1.0, 1.0, 0.0, 0.0},
        {-1.0, 0.0, 1.0, 0.0},
        {-1.0, 0.0, 0.0, 1.0}
    }}
    , paramBeta_{{
        {0.0, 0.25, 0.0,  0.0},   // stage 1: a_11 = 1/4
        {0.0, 0.5,  0.25, 0.0},   // stage 2: a_21 = 1/2, a_22 = 1/4
        {0.0, 0.5,  0.5,  0.0}    // update:  b_1  = 1/2, b_2  = 1/2
    }}
    , paramD_{{0.0, 0.25, 0.75, 1.0}}
    {}

    bool implicit () const final
    { return true; }

    std::size_t numStages () const final
    { return 3; }

    Scalar temporalWeight (std::size_t i, std::size_t k) const final
    { return paramAlpha_[i-1][k]; }

    Scalar spatialWeight (std::size_t i, std::size_t k) const final
    { return paramBeta_[i-1][k]; }

    Scalar timeStepWeight (std::size_t k) const final
    { return paramD_[k]; }

    std::string name () const final
    { return "symplectic diagonally implicit Runge-Kutta 2nd order (Qin-Zhang)"; }

private:
    std::array<std::array<Scalar, 4>, 3> paramAlpha_;
    std::array<std::array<Scalar, 4>, 3> paramBeta_;
    std::array<Scalar, 4> paramD_;
};

} // end namespace MultiStage
} // end namespace Dumux::Experimental

#endif
