// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Additional source term of the radial momentum balance in rotationally symmetric problems
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_AXISYMMETRIC_SOURCE_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_AXISYMMETRIC_SOURCE_HH

namespace Dumux::Detail {

/*!
 * \ingroup NavierStokesModel
 * \brief Source density of the radial momentum balance of two-dimensional problems with
 *        rotational extrusion (axisymmetric flow without swirl)
 *
 * With the viscous fluxes and the pressure integrated over the faces of the rotated control
 * volumes (as \f$ \boldsymbol{\tau}\mathbf{n} \f$ and \f$ p\mathbf{n} \f$), the radial momentum
 * balance in cylindrical coordinates contains the additional source density
 * \f$ -\tau_{\theta\theta}/r + p/r \f$, see Ferziger/Peric: Computational methods for Fluid
 * Dynamics (2020), https://doi.org/10.1007/978-3-319-99693-6, Chapter 9.9 and Eq. (9.81). The hoop
 * stress is \f$ \tau_{\theta\theta} = 2\mu u_r/r \f$ for the viscous stress
 * \f$ \mu(\nabla\mathbf{u} + \nabla\mathbf{u}^T) \f$ and \f$ \tau_{\theta\theta} = \mu u_r/r \f$
 * for \f$ \mu\nabla\mathbf{u} \f$.
 *
 * \param radius The distance \f$ r \f$ to the rotation axis
 * \param radialVelocity The radial velocity \f$ u_r \f$
 * \param viscosity The dynamic viscosity \f$ \mu \f$
 * \param pressure The pressure \f$ p \f$
 * \param unsymmetrizedVelocityGradient If the viscous stress is \f$ \mu\nabla\mathbf{u} \f$
 */
template<class Scalar>
Scalar axisymmetricRadialMomentumSource(Scalar radius,
                                        Scalar radialVelocity,
                                        Scalar viscosity,
                                        Scalar pressure,
                                        bool unsymmetrizedVelocityGradient)
{
    const Scalar hoopStressFactor = unsymmetrizedVelocityGradient ? 1.0 : 2.0;
    return -hoopStressFactor*viscosity*radialVelocity/(radius*radius) + pressure/radius;
}

} // end namespace Dumux::Detail

#endif
