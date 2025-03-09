// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Elastic
 * \brief Deformation of a solid body using the theory of linear elasticity (small deformations)
 *
 * This model describes the deformation of a solid body using the theory of linear elasticity.
 * The momentum balance equation of the solid body can be expressed by
 \f[
 \nabla\cdot\boldsymbol{\sigma} + \rho \mathbf{g} + \mathbf{f} = \rho\ddot{\mathbf{u}},
 \f]
 * where \f$ \boldsymbol{\sigma} \f$ is the stress tensor, \f$ \rho \f$ is the density of
 * the solid, \f$ \mathbf{f} \f$ in \f$ \mathrm{N/m^3} \f$  is the external force acting on the body per unit volume (e.g. magnetism),
 * and \f$ \mathbf{u} = \mathbf{x} - \mathbf{x}_{\mathrm{initial}} \f$ is the displacement,
 * defined as the difference in material points \f$ \mathbf{x} \f$ and \f$ \mathbf{x}_{\mathrm{initial}} \f$
 * in the deformed and undeformed (initial) state, respectively. The model assumes quasi-static conditions,
 * that is, the above momentum balance equation is solved under the assumption that the acceleration term
 * \f$ \rho\ddot{\mathbf{u}} \approx 0\f$.
 *
 * Per default, Hookes' Law is used for expressing the stress tensor \f$ \boldsymbol{\sigma} \f$ as a function of the
 * displacement:
 \f[
 \boldsymbol{\sigma} = \lambda\mathrm{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2G \boldsymbol{\varepsilon},
 \f]
 * with
 \f[
 \boldsymbol{\varepsilon} = \frac{1}{2} \left[ \nabla\mathbf{u} + (\nabla\mathbf{u})^{\mathrm{T}} \right].
 \f]
 *
 * Primary variables are the displacements in each direction \f$ \mathbf{u} \f$.
 * Gravity can be enabled or disabled via a runtime parameter.
 *
 */
#ifndef DUMUX_GEOMECHANICS_ELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_ELASTIC_MODEL_HH

#include <dumux/solidmechanics/elastic/model.hh>
#warning "This header is deprecated and will be removed after 3.10. Use dumux/solidmechanics/elastic/model.hh."

 #endif
