// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Parameters of linear multi-step time integration methods
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTEP_METHODS_HH
#define DUMUX_TIMESTEPPING_MULTISTEP_METHODS_HH

#include <array>
#include <cmath>
#include <span>
#include <string>
#include <algorithm>

#include <dune/common/exceptions.hh>

#include <dumux/experimental/timestepping/timeintegrationmethod.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \brief Abstract interface for linear multi-step method parameters
 *
 * We consider the semi-discrete problem
 *
 * \f[
 *   \frac{\partial M(x, t)}{\partial t} + R(x, t) = 0
 * \f]
 *
 * (see ImplicitSolveParams). A \f$ k \f$-step method computes the solution
 * \f$ x^{n+1} \f$ at \f$ t^{n+1} = t^n + \Delta t^n \f$ from the solutions at
 * the \f$ k \f$ previous time levels by solving
 *
 * \f[
 *   \sum_{j=0}^{k} \left[ \alpha_j M\left(x^{n+1-j}, t^{n+1-j}\right)
 *     + \beta_j \Delta t^n R\left(x^{n+1-j}, t^{n+1-j}\right) \right] = 0.
 * \f]
 *
 * For variable time step sizes the weights \f$ \alpha_j, \beta_j \f$ depend on the
 * step sizes \f$ \Delta t^n, \Delta t^{n-1}, \ldots, \Delta t^{n+1-k} \f$.
 * The method also decides how many previous time levels to use for given step
 * sizes, which allows for a lower order at the start of a simulation and for
 * step size sequences outside the stability limits of the full order.
 */
template<class Scalar>
class MultiStepMethod : public TimeIntegrationMethod<Scalar>
{
public:
    //! The largest number of previous time levels \f$ k \f$ the method uses
    virtual std::size_t numSteps () const = 0;

    //! Whether the spatial operator at previous time levels enters (\f$ \beta_j \neq 0 \f$ for some \f$ j \geq 1 \f$)
    virtual bool spatialHistory () const = 0;

    /*!
     * \brief The number of previous time levels to use for the next step
     * \param stepSizes The size \f$ \Delta t^n \f$ of the next step followed by the sizes
     *        \f$ \Delta t^{n-1}, \Delta t^{n-2}, \ldots \f$ of the previous steps, one for each
     *        previous time level that is available in addition to \f$ t^n \f$
     * \return a number \f$ k \f$ with \f$ 1 \leq k \leq \f$ `stepSizes.size()`
     */
    virtual std::size_t numStepsToUse (std::span<const Scalar> stepSizes) const = 0;

    /*!
     * \brief The weights for the given step sizes
     * \param stepSizes The step sizes \f$ \Delta t^n, \ldots, \Delta t^{n+1-k} \f$ (size \f$ k \f$)
     * \param temporal The temporal weights \f$ \alpha_0, \ldots, \alpha_k \f$ (output, size \f$ k+1 \f$)
     * \param spatial The spatial weights \f$ \beta_0, \ldots, \beta_k \f$ (output, size \f$ k+1 \f$)
     */
    virtual void weights (std::span<const Scalar> stepSizes,
                          std::span<Scalar> temporal,
                          std::span<Scalar> spatial) const = 0;
};

//! Multi-step time stepping scheme implementations
namespace MultiStep {

/*!
 * \brief Backward differentiation formula (BDF) of order \f$ k \leq 5 \f$
 *
 * The storage term is differentiated by the interpolation polynomial through the
 * time levels \f$ t^{n+1}, \ldots, t^{n+1-k} \f$ and the spatial operator is evaluated
 * at the new time level only (\f$ \beta = (1, 0, \ldots, 0) \f$):
 *
 * \f[
 *   \alpha_j = \Delta t^n \, \ell_j'(t^{n+1}), \quad j = 0, \ldots, k,
 * \f]
 *
 * where \f$ \ell_j \f$ are the Lagrange basis polynomials of the time levels.
 * Because the spatial operator is not evaluated at previous time levels, errors in
 * equations without storage term (constraints) do not propagate between steps.
 *
 * BDF1 (implicit Euler) and BDF2 are A-stable. BDF2 with variable step sizes is
 * zero-stable for step size ratios \f$ \omega = \Delta t^n/\Delta t^{n-1} < 1 + \sqrt{2} \f$
 * \cite Grigorieff1983. BDF3 to BDF5 are \f$ A(\alpha) \f$-stable with
 * \f$ \alpha = 86.03^\circ, 73.35^\circ, 51.84^\circ \f$ \cite HairerWanner1996 (Section V.2)
 * and are only used for equal step sizes.
 *
 * The number of time levels used for a step is the largest \f$ k \f$ not larger than the
 * order and the number of available levels for which the last \f$ k \f$ step sizes are
 * equal (for \f$ k \geq 3 \f$) or the step size ratio does not exceed the given maximum
 * (for \f$ k = 2 \f$). Otherwise the step is taken with BDF1.
 */
template<class Scalar>
class BDF final : public MultiStepMethod<Scalar>
{
public:
    /*!
     * \brief Construct a BDF method
     * \param order The order \f$ k \f$ (between 1 and 5)
     * \param maxStepSizeRatio The largest ratio \f$ \Delta t^n/\Delta t^{n-1} \f$ for which BDF2 is used
     * \param stepSizeTolerance The relative tolerance for considering step sizes equal
     */
    explicit BDF(const std::size_t order,
                 const Scalar maxStepSizeRatio = 2.0,
                 const Scalar stepSizeTolerance = 1e-8)
    : order_(order)
    , maxStepSizeRatio_(maxStepSizeRatio)
    , stepSizeTolerance_(stepSizeTolerance)
    {
        if (order_ < 1 || order_ > 5)
            DUNE_THROW(Dune::InvalidStateException, "BDF methods are implemented for orders 1 to 5, not " << order_);

        using std::sqrt;
        if (!(maxStepSizeRatio_ > 0.0 && maxStepSizeRatio_ < 1.0 + sqrt(2.0)))
            DUNE_THROW(Dune::InvalidStateException,
                "The maximum step size ratio of BDF2 has to be in (0, 1 + sqrt(2)), not " << maxStepSizeRatio_);
    }

    bool implicit () const final
    { return true; }

    std::size_t numSteps () const final
    { return order_; }

    bool spatialHistory () const final
    { return false; }

    std::size_t numStepsToUse (std::span<const Scalar> stepSizes) const final
    {
        if (stepSizes.empty())
            DUNE_THROW(Dune::InvalidStateException, "At least the size of the next step is required");

        auto k = std::min(order_, stepSizes.size());
        while (k > 2 && !equalStepSizes_(stepSizes.first(k)))
            --k;

        if (k == 2 && stepSizes[0] > maxStepSizeRatio_*stepSizes[1])
            k = 1;

        return k;
    }

    void weights (std::span<const Scalar> stepSizes,
                  std::span<Scalar> temporal,
                  std::span<Scalar> spatial) const final
    {
        const auto k = stepSizes.size();
        if (k < 1 || k > order_ || temporal.size() != k+1 || spatial.size() != k+1)
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of step sizes or weights");

        // time levels relative to the new one, tau_0 = 0 > tau_1 > ... > tau_k
        std::array<Scalar, 6> tau{};
        for (std::size_t j = 1; j <= k; ++j)
            tau[j] = tau[j-1] - stepSizes[j-1];

        const auto dt = stepSizes[0];
        temporal[0] = 0.0;
        for (std::size_t m = 1; m <= k; ++m)
            temporal[0] -= dt/tau[m];

        for (std::size_t j = 1; j <= k; ++j)
        {
            Scalar numerator = 1.0, denominator = 1.0;
            for (std::size_t m = 0; m <= k; ++m)
            {
                if (m == j)
                    continue;
                if (m != 0)
                    numerator *= -tau[m];
                denominator *= tau[j] - tau[m];
            }
            temporal[j] = dt*numerator/denominator;
        }

        std::fill(spatial.begin(), spatial.end(), 0.0);
        spatial[0] = 1.0;
    }

    std::string name () const final
    { return "BDF" + std::to_string(order_); }

private:
    bool equalStepSizes_(std::span<const Scalar> stepSizes) const
    {
        using std::abs;
        return std::all_of(stepSizes.begin(), stepSizes.end(), [&](const Scalar dt)
        { return abs(dt - stepSizes[0]) <= stepSizeTolerance_*stepSizes[0]; });
    }

    std::size_t order_;
    Scalar maxStepSizeRatio_;
    Scalar stepSizeTolerance_;
};

} // end namespace MultiStep
} // end namespace Dumux::Experimental

#endif
