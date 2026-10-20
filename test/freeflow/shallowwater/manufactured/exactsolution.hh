// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief A smooth manufactured solution of the shallow water equations with a non-flat bed
 */
#ifndef DUMUX_SHALLOWWATER_MANUFACTURED_EXACT_SOLUTION_HH
#define DUMUX_SHALLOWWATER_MANUFACTURED_EXACT_SOLUTION_HH

#include <array>
#include <cmath>

#include <dune/common/fvector.hh>

namespace Dumux::ShallowWaterTest {

/*!
 * \ingroup ShallowWaterTests
 * \brief A smooth manufactured solution of the shallow water equations with a non-flat bed
 *
 * On the unit square, with the frequency \f$ \omega \f$,
 * \f[
 *   h = 1 + \frac{1}{10} \sin(\pi x + \omega t) \cos(\pi y), \quad
 *   u = \frac{1}{2} + \frac{1}{10} \cos(\pi x) \sin(2\pi y + \omega t), \quad
 *   v = \frac{1}{5} \sin(\pi x) \sin(\pi y) \cos(\omega t), \quad
 *   z = \frac{1}{5} \exp\left(-10\left((x - \tfrac{1}{2})^2 + (y - \tfrac{1}{2})^2\right)\right),
 * \f]
 * a subcritical flow (Froude number below 0.21 for \f$ g = 9.81 \f$), steady for \f$ \omega = 0 \f$.
 * The sources are the residuals of the equations, \f$ \partial_t h + \nabla\cdot(h\mathbf{u}) \f$ for the mass
 * and \f$ \partial_t (h\mathbf{u}) + \nabla\cdot(h\mathbf{u}\otimes\mathbf{u}) + \nabla(g h^2/2) + g h \nabla z \f$
 * for the momentum.
 */
template<class Scalar>
class ManufacturedSolution
{
    using GlobalPosition = Dune::FieldVector<Scalar, 2>;

public:
    //! water depth, velocity components and bed surface at a position and time for the frequency omega
    static std::array<Scalar, 4> state(const GlobalPosition& pos, const Scalar t = 0.0, const Scalar omega = 0.0)
    {
        using std::sin; using std::cos; using std::exp;
        const Scalar x = pos[0], y = pos[1];
        return {{ 1.0 + 0.1*sin(M_PI*x + omega*t)*cos(M_PI*y),
                  0.5 + 0.1*cos(M_PI*x)*sin(2.0*M_PI*y + omega*t),
                  0.2*sin(M_PI*x)*sin(M_PI*y)*cos(omega*t),
                  0.2*exp(-10.0*((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5))) }};
    }

    //! the sources of mass and momentum (per unit area) for the gravitational acceleration g at a time
    //! for the frequency omega
    static std::array<Scalar, 3> source(const GlobalPosition& pos, const Scalar g,
                                        const Scalar t = 0.0, const Scalar omega = 0.0)
    {
        using std::sin; using std::cos;
        const Scalar x = pos[0], y = pos[1];
        const Scalar phaseH = M_PI*x + omega*t, phaseU = 2.0*M_PI*y + omega*t;
        const auto [h, u, v, z] = state(pos, t, omega);
        const Scalar ht = 0.1*omega*cos(phaseH)*cos(M_PI*y);
        const Scalar hx = 0.1*M_PI*cos(phaseH)*cos(M_PI*y), hy = -0.1*M_PI*sin(phaseH)*sin(M_PI*y);
        const Scalar ut = 0.1*omega*cos(M_PI*x)*cos(phaseU);
        const Scalar ux = -0.1*M_PI*sin(M_PI*x)*sin(phaseU), uy = 0.2*M_PI*cos(M_PI*x)*cos(phaseU);
        const Scalar vt = -0.2*omega*sin(M_PI*x)*sin(M_PI*y)*sin(omega*t);
        const Scalar vx = 0.2*M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(omega*t), vy = 0.2*M_PI*sin(M_PI*x)*cos(M_PI*y)*cos(omega*t);
        const Scalar zx = -20.0*(x - 0.5)*z, zy = -20.0*(y - 0.5)*z;

        return {{ ht + hx*u + h*ux + hy*v + h*vy,
                  ht*u + h*ut + hx*u*u + 2.0*h*u*ux + g*h*hx + hy*u*v + h*uy*v + h*u*vy + g*h*zx,
                  ht*v + h*vt + hx*u*v + h*ux*v + h*u*vx + hy*v*v + 2.0*h*v*vy + g*h*hy + g*h*zy }};
    }
};

} // end namespace Dumux::ShallowWaterTest

#endif
