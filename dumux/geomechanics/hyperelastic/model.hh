// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Hyperelastic model
 * \brief Deformation of a solid body using the theory of (nonlinear) elasticity (large deformations)
 *
 * This model describes the deformation of a solid body using the finite strain theory.
 * The displacement vector of a material point following the deformation function
 * \f$ \mathbf{x} = \chi(\mathbf{X}) \f$ is given by the difference of current position
 * and position in the reference configuration:
 \f[
 \mathbf{d} = \mathbf{x} - \mathbf{X}.
 \f]
 * An important kinematic quantity characterizing the deformation is the deformation gradient
 \f[
 \mathbf{F} = \frac{\partial \mathbf{x}}{\partial \mathbf{X}} = \mathbf{I} + \nabla \mathbf{d}.
 \f]
 * The equilibrium equation (nonlinear elastostatics) in the reference frame is given by
 \f[
 - \nabla\cdot\mathbf{P} = \mathbf{f}
 \f]
 * where \f$ \mathbf{f} \f$ in \f$ \mathrm{N/m^3} \f$ is the external force acting
 * on the body per unit volume and \f$ \mathbf{P} = \mathbf{P}(\mathbf{F}) \f$ is
 * the first Piola-Kirchhoff stress tensor which relates the stress in the current
 * and the reference configuration. It is related to the Cauchy stress,
 * \f$ \mathbf{P} = J \boldsymbol{\sigma} \mathbf{F}^{-T} \f$,
 * where \f$ J = \operatorname{det}{\mathbf{F}} \f$.
 *
 * A suitable constitutive law for \f$ \mathbf{P} = \mathbf{P}(\mathbf{F}) \f$ completes
 * the model. All hyperelastic materials can be described in terms of a strain energy
 * density function \f$ \psi(\mathbf{F}) \f$ @cite ogden1997
 * and the first Piola-Kirchhoff stress tensor can be computed as
 \f[
 \mathbf{P} = \frac{\partial \psi}{\partial \mathbf{F}}
            = \mathbf{F}2\frac{\partial \psi}{\partial \mathbf{C}},
 \f]
 * where \f$ \mathbf{C} = \mathbf{F}^T\mathbf{F} \f$ is the right Cauchy-Green tensor
 * and the term \f$ \mathbf{S} = 2 (\partial \psi / \partial \mathbf{C}) \f$ is the
 * second Piola-Kirchhoff stress tensor.
 * This model expects the user problem implementation to provide a function
 * `firstPiolaKirchhoffStressTensor(F)` implementing the constitutive law.
 */
#ifndef DUMUX_GEOMECHANICS_HYPERELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_HYPERELASTIC_MODEL_HH

#include <dumux/solidmechanics/hyperelastic/model.hh>
#warning "This header is deprecated and will be removed after 3.10. Use dumux/solidmechanics/hyperelastic/model.hh."

#endif
