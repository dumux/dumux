// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Analytical solutions and rhs for the different test cases
 */

#ifndef DUMUX_CONVERGENCE_TEST_ANALYTICAL_SOLUTIONS_HH
#define DUMUX_CONVERGENCE_TEST_ANALYTICAL_SOLUTIONS_HH

#include <cmath>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux::Solution::DarcyStokes {

/////////////////////////////////////////////////////////////////////////////////////////////////
// see Rybak et al., 2015: "Multirate time integration for coupled
//   saturated/unsaturated porous medium and free flow systems"
/////////////////////////////////////////////////////////////////////////////////////////////////
namespace Rybak {

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> darcy(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];

    using std::exp; using std::sin; using std::cos;
    sol[0] = -0.5*M_PI*y*y*cos(M_PI*x);
    sol[1] = -1.0*y*sin(M_PI*x);
    sol[2] = 0.5*y*y*sin(M_PI*x);
    return sol;
}

// 0: mass
Dune::FieldVector<double, 1> darcyRHS(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];
    using std::sin;
    return { (0.5*M_PI*y*M_PI*y - 1)*sin(M_PI*x) };
}

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> stokes(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];

    using std::sin; using std::cos;
    sol[0] = -cos(M_PI*x)*sin(M_PI*y);
    sol[1] = sin(M_PI*x)*cos(M_PI*y);
    sol[2] = 0.5*y*sin(M_PI*x);
    return sol;
}

// 0: mom-x | 1: mom-y | 2: mass
Dune::FieldVector<double, 3> stokesRHS(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];
    Dune::FieldVector<double, 3> source(0.0);
    using std::sin; using std::cos;
    source[0] = 0.5*M_PI*y*cos(M_PI*x) - 2*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y);
    source[1] = 0.5*sin(M_PI*x) + 2*M_PI*M_PI*sin(M_PI*x)*cos(M_PI*y);
    return source;
}

// sigma_ff = -sym(grad(v)) + pI
Dune::FieldMatrix<double, 2, 2> stokesStress(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];

    using std::exp; using std::sin; using std::cos;

    Dune::FieldMatrix<double, 2, 2> stress(0.0);
    stress[0][0] = (0.5*y - 2*M_PI*sin(M_PI*y))*sin(M_PI*x);
    stress[1][0] = 0.0;
    stress[0][1] = stress[1][0]; // symmetric
    stress[1][1] = (0.5*y - 2*M_PI*sin(M_PI*y))*sin(M_PI*x);
    return stress;
}

} // end namespace Rybak

/////////////////////////////////////////////////////////////////////////////////////////////////
// see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
// Example 1
/////////////////////////////////////////////////////////////////////////////////////////////////
namespace ShiueOne {

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> darcy(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];

    using std::exp; using std::sin; using std::cos;
    sol[0] = M_PI*(-exp(1)*y + exp(y))*sin(M_PI*x);
    sol[1] = (exp(1) - exp(y))*cos(M_PI*x);
    sol[2] = (exp(y) - y*exp(1)) * cos(M_PI*x);

    return sol;
}

// 0: mass
Dune::FieldVector<double, 1> darcyRHS(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];
    using std::exp; using std::sin; using std::cos;
    return { cos(M_PI*x) * (M_PI*M_PI*(exp(y) - y*exp(1)) - exp(y)) };
}

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> stokes(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];

    using std::exp; using std::sin; using std::cos;
    sol[0] = -1/M_PI * exp(y) * sin(M_PI*x);
    sol[1] = (exp(y) - exp(1)) * cos(M_PI*x);
    sol[2] = 2*exp(y) * cos(M_PI*x);
    return sol;
}

// 0: mom-x | 1: mom-y | 2: mass
Dune::FieldVector<double, 3> stokesRHS(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];
    using std::exp; using std::sin; using std::cos;
    Dune::FieldVector<double, 3> source(0.0);
    source[0] = exp(y)*sin(M_PI*x) * (1/M_PI -3*M_PI);
    source[1] = cos(M_PI*x) * (M_PI*M_PI*(exp(y)- exp(1)) + exp(y));
    return source;
}

// sigma_ff = -sym(grad(v)) + pI
Dune::FieldMatrix<double, 2, 2> stokesStress(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];

    using std::exp; using std::sin; using std::cos;

    Dune::FieldMatrix<double, 2, 2> stress(0.0);
    stress[0][0] = 4.0*exp(y)*cos(M_PI*x);
    stress[1][0] = (M_PI*M_PI*(exp(y) - exp(1)) + 1.0*exp(y))*sin(M_PI*x)/M_PI;
    stress[0][1] = stress[1][0]; // symmetric
    stress[1][1] = 0.0;
    return stress;
}

} // end namespace ShiueExOne

/////////////////////////////////////////////////////////////////////////////////////////////////
// see Shiue et al., 2018: "Convergence of the MAC Scheme for the Stokes/Darcy Coupling Problem"
// Example 2
/////////////////////////////////////////////////////////////////////////////////////////////////
namespace ShiueTwo {

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> darcy(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];

    sol[0] = x*(y - 1.0) + (x - 1.0)*(y - 1.0) - 2.0;
    sol[1] = x*(x - 1.0) - 1.0*(y - 1.0)*(y - 1.0) - 2.0;
    sol[2] = x*(1.0-x)*(y-1.0) + (y-1.0)*(y-1.0)*(y-1.0)/3.0 + 2.0*x + 2.0*y + 4.0;

    return sol;
}

// 0: mass
Dune::FieldVector<double, 1> darcyRHS(const Dune::FieldVector<double, 2>& globalPos)
{ return { 0.0 }; }

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> stokes(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];

    sol[0] = (y-1.0)*(y-1.0) + x*(y-1.0) + 3.0*x - 1.0;
    sol[1] = x*(x-1.0) - 0.5*(y-1.0)*(y-1.0) - 3.0*y + 1.0;
    sol[2] = 2.0*x + y - 1.0;
    return sol;
}

// 0: mom-x | 1: mom-y | 2: mass
Dune::FieldVector<double, 3> stokesRHS(const Dune::FieldVector<double, 2>& globalPos)
{ return { 0.0, 0.0, 0.0 }; }

// sigma_ff = -sym(grad(v)) + pI
Dune::FieldMatrix<double, 2, 2> stokesStress(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];

    Dune::FieldMatrix<double, 2, 2> stress(0.0);
    stress[0][0] = 2.0*x - y - 5.0;
    stress[1][0] = -3*x - 2*y + 3.0;
    stress[0][1] = stress[1][0]; // symmetric
    stress[1][1] = 2.0*x + 3.0*y + 3.0;
    return stress;
}

} // end namespace ShiueTwo

/////////////////////////////////////////////////////////////////////////////////////////////////
// see Schneider et al., 2019: "Coupling staggered-grid and MPFA finite volume methods for
//   free flow/porous-medium flow problems (with c = 0)"
/////////////////////////////////////////////////////////////////////////////////////////////////
namespace Schneider {

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> darcy(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];
    static constexpr double omega = M_PI;
    static constexpr double c = 0.0;
    using std::exp; using std::sin; using std::cos;
    const double sinOmegaX = sin(omega*x);
    const double cosOmegaX = cos(omega*x);
    static const double expTwo = exp(2);
    const double expYPlusOne = exp(y+1);

    sol[0] = c/(2*omega)*expYPlusOne*sinOmegaX*sinOmegaX
              -omega*(expYPlusOne + 2 - expTwo)*cosOmegaX;
    sol[1] = (0.5*c*(expYPlusOne + 2 - expTwo)*cosOmegaX
              -(c*cosOmegaX + 1)*exp(y-1))*sinOmegaX;
    sol[2] = (expYPlusOne + 2 - expTwo)*sinOmegaX;

    return sol;
}

// 0: mass
Dune::FieldVector<double, 1> darcyRHS(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];
    using std::exp; using std::sin; using std::cos;
    static constexpr double omega = M_PI;
    static constexpr double c = 0.0;
    const double cosOmegaX = cos(omega*x);
    static const double expTwo = exp(2);
    const double expYPlusOne = exp(y+1);

    return { sin(omega*x)*(-(c*cosOmegaX + 1)*exp(y - 1)
            + 1.5*c*expYPlusOne*cosOmegaX
            + omega*omega*(expYPlusOne - expTwo + 2)) };
}

// 0: velocity-x | 1: velocity-y | 2: pressure
Dune::FieldVector<double, 3> stokes(const Dune::FieldVector<double, 2>& globalPos)
{
    Dune::FieldVector<double, 3> sol(0.0);
    const double x = globalPos[0];
    const double y = globalPos[1];
    using std::sin;
    static constexpr double omega = M_PI;
    const double sinOmegaX = sin(omega*x);

    sol[0] = y;
    sol[1] = -y*sinOmegaX;
    sol[2] = -y*y*sinOmegaX*sinOmegaX;
    return sol;
}

// 0: mom-x | 1: mom-y | 2: mass
Dune::FieldVector<double, 3> stokesRHS(const Dune::FieldVector<double, 2>& globalPos)
{
    const double x = globalPos[0];
    const double y = globalPos[1];
    using std::exp; using std::sin; using std::cos;
    static constexpr double omega = M_PI;
    const double sinOmegaX = sin(omega*x);
    const double cosOmegaX = cos(omega*x);

    Dune::FieldVector<double, 3> source(0.0);
    source[0] = -2*omega*y*y*sinOmegaX*cosOmegaX
                                            -2*y*sinOmegaX + omega*cosOmegaX;
    source[1] = -omega*y*y*cosOmegaX - omega*omega*y*sinOmegaX;
    source[2] = -sinOmegaX;
    return source;
}

// sigma_ff = vvT + sym(grad(v)) + pI
Dune::FieldMatrix<double, 2, 2> stokesStress(const Dune::FieldVector<double, 2>& globalPos)
{
    using std::sin; using std::cos;
    const double x = globalPos[0];
    const double y = globalPos[1];
    static constexpr double omega = M_PI;
    const double sinOmegaX = sin(omega*x);
    const double cosOmegaX = cos(omega*x);

    Dune::FieldMatrix<double, 2, 2> stress(0.0);
    stress[0][0] = y*y * cosOmegaX*cosOmegaX;
    stress[0][1] = -y*y*sinOmegaX + omega*y*cosOmegaX - 1.0;
    stress[1][0] = stress[0][1]; // symmetric
    stress[1][1] = 2*sinOmegaX;
    return stress;
}

} // end namespace Schneider

} // end namespace Dumux::Solution::DarcyStokes

#endif
