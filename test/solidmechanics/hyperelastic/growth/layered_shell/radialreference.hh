// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Radially symmetric solution of a two-layer hollow cylinder (plane strain) or sphere with
 *        isotropic growth in each layer
 *
 * The material point at radius \f$ R \f$ moves to radius \f$ r(R) \f$, with the radial stretch
 * \f$ \lambda_r = r' \f$ and the \f$ n = d - 1 \f$ hoop stretches \f$ \lambda_\theta = r/R \f$.
 * In a layer with growth stretch \f$ \vartheta \f$, the principal first Piola-Kirchhoff stresses
 * of the compressible neo-Hookean law
 * \f$ \mathbf{P}_\mathrm{e} = \mu \mathbf{F}_\mathrm{e} + (\lambda \ln J_\mathrm{e} - \mu)\mathbf{F}_\mathrm{e}^{-T} \f$
 * pulled back with \f$ J_\mathrm{g} \mathbf{F}_\mathrm{g}^{-T} = \vartheta^{n} \mathbf{I} \f$ are
 * \f[ P_i = \vartheta^n \left( \mu a_i + \frac{\lambda \ln J_\mathrm{e} - \mu}{a_i} \right), \quad
 *     a_i = \lambda_i / \vartheta, \quad J_\mathrm{e} = a_r a_\theta^n. \f]
 * Radial equilibrium \f$ P_r' + n (P_r - P_\theta)/R = 0 \f$ is integrated for the state
 * \f$ (r, P_r) \f$, which is continuous across the layer interface, with the classical Runge-Kutta
 * method. The unknown \f$ r(R_i) \f$ is found by the secant method such that both surfaces are
 * traction free.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_RADIAL_REFERENCE_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_RADIAL_REFERENCE_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <dune/common/exceptions.hh>

namespace Dumux {

template<class Scalar>
class LayeredShellRadialReference
{
public:
    struct Layer { Scalar growthStretch, shearModulus, lameLambda; };

    LayeredShellRadialReference(int numHoopDirections,
                                std::array<Scalar, 3> radii,
                                std::array<Layer, 2> layers,
                                int stepsPerLayer = 20000)
    : n_(numHoopDirections), radii_(radii), layers_(layers), steps_(stepsPerLayer)
    {
        Scalar r0 = layers_[0].growthStretch*radii_[0];
        Scalar r1 = 1.01*r0;
        Scalar f0 = integrate_(r0, false);
        for (int iter = 0; iter < 50; ++iter)
        {
            const Scalar f1 = integrate_(r1, false);
            if (std::abs(f1) < 1e-13)
                break;
            const Scalar r2 = r1 - f1*(r1 - r0)/(f1 - f0);
            r0 = r1; f0 = f1; r1 = r2;
        }

        const Scalar residual = integrate_(r1, true);
        if (std::abs(residual) > 1e-10)
            DUNE_THROW(Dune::MathError, "Radial reference: outer surface not traction free, P_r = " << residual);
    }

    //! The current radius \f$ r(R) \f$ (cubic Hermite interpolation with \f$ r' = \lambda_r \f$)
    Scalar radius(const Scalar R) const
    {
        const auto [k, t, h] = locate_(R);
        const Scalar t2 = t*t, t3 = t2*t;
        return (2*t3 - 3*t2 + 1)*r_[k] + (t3 - 2*t2 + t)*h*stretch_[k]
               + (-2*t3 + 3*t2)*r_[k+1] + (t3 - t2)*h*stretch_[k+1];
    }

    //! The radial stretch \f$ \lambda_r = r'(R) \f$ (linear interpolation)
    Scalar radialStretch(const Scalar R) const
    {
        const auto [k, t, h] = locate_(R);
        return (1.0 - t)*stretch_[k] + t*stretch_[k+1];
    }

private:
    Scalar radialStress_(Scalar lr, Scalar lt, const Layer& l) const
    {
        const Scalar a = lr/l.growthStretch, b = lt/l.growthStretch;
        const Scalar lnJe = std::log(a*std::pow(b, n_));
        return std::pow(l.growthStretch, n_)*(l.shearModulus*a + (l.lameLambda*lnJe - l.shearModulus)/a);
    }

    Scalar hoopStress_(Scalar lr, Scalar lt, const Layer& l) const
    {
        const Scalar a = lr/l.growthStretch, b = lt/l.growthStretch;
        const Scalar lnJe = std::log(a*std::pow(b, n_));
        return std::pow(l.growthStretch, n_)*(l.shearModulus*b + (l.lameLambda*lnJe - l.shearModulus)/b);
    }

    // the radial stretch for given radial stress and hoop stretch (Newton's method)
    Scalar radialStretch_(Scalar Pr, Scalar lt, const Layer& l, Scalar guess) const
    {
        Scalar lr = guess;
        for (int iter = 0; iter < 50; ++iter)
        {
            const Scalar a = lr/l.growthStretch, b = lt/l.growthStretch;
            const Scalar lnJe = std::log(a*std::pow(b, n_));
            const Scalar f = radialStress_(lr, lt, l) - Pr;
            const Scalar df = std::pow(l.growthStretch, n_ - 1)
                              *(l.shearModulus + (l.lameLambda*(1.0 - lnJe) + l.shearModulus)/(a*a));
            lr -= f/df;
            if (std::abs(f) < 1e-15)
                break;
        }
        return lr;
    }

    // right-hand side of (r, P_r)' and the radial stretch
    std::array<Scalar, 3> rhs_(Scalar R, Scalar r, Scalar Pr, const Layer& l, Scalar guess) const
    {
        const Scalar lt = r/R;
        const Scalar lr = radialStretch_(Pr, lt, l, guess);
        return {lr, -n_*(Pr - hoopStress_(lr, lt, l))/R, lr};
    }

    // integrates from the traction-free inner surface and returns P_r at the outer surface
    Scalar integrate_(Scalar rInner, bool store)
    {
        if (store)
        {
            R_.clear(); r_.clear(); stretch_.clear();
        }

        Scalar r = rInner, Pr = 0.0;
        Scalar guess = layers_[0].growthStretch;
        for (int layerIdx = 0; layerIdx < 2; ++layerIdx)
        {
            const auto& l = layers_[layerIdx];
            const Scalar Ra = radii_[layerIdx], Rb = radii_[layerIdx+1];
            const int steps = store ? steps_ : steps_/10;
            const Scalar h = (Rb - Ra)/steps;
            for (int step = 0; step <= steps; ++step)
            {
                const Scalar R = Ra + step*h;
                const auto k1 = rhs_(R, r, Pr, l, guess);
                guess = k1[2];
                if (store)
                {
                    R_.push_back(R); r_.push_back(r); stretch_.push_back(k1[2]);
                }
                if (step == steps)
                    break;

                const auto k2 = rhs_(R + 0.5*h, r + 0.5*h*k1[0], Pr + 0.5*h*k1[1], l, guess);
                const auto k3 = rhs_(R + 0.5*h, r + 0.5*h*k2[0], Pr + 0.5*h*k2[1], l, guess);
                const auto k4 = rhs_(R + h, r + h*k3[0], Pr + h*k3[1], l, guess);
                r += h/6.0*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
                Pr += h/6.0*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
            }
        }
        return Pr;
    }

    // interval index, local coordinate and length; the layers are stored one after the other
    struct Location { std::size_t k; Scalar t; Scalar h; };
    Location locate_(Scalar R) const
    {
        R = std::clamp(R, R_.front(), R_.back());
        const std::size_t offset = (R < radii_[1]) ? 0 : steps_ + 1;
        const Scalar Ra = radii_[R < radii_[1] ? 0 : 1];
        const Scalar Rb = radii_[R < radii_[1] ? 1 : 2];
        const Scalar h = (Rb - Ra)/steps_;
        const std::size_t local = std::min<std::size_t>(static_cast<std::size_t>((R - Ra)/h), steps_ - 1);
        return {offset + local, (R - R_[offset + local])/h, h};
    }

    int n_;
    std::array<Scalar, 3> radii_;
    std::array<Layer, 2> layers_;
    int steps_;
    std::vector<Scalar> R_, r_, stretch_;
};

} // end namespace Dumux

#endif
