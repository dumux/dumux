// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElastic
 * \brief A poroelastic geomechanical model
 *
 * This model describes the deformation of a porous medium using the theory of linear poroelasticity.
 * The momentum balance equation of a porous medium can be expressed by
 \f[
 \nabla\cdot\boldsymbol{\sigma_{\mathrm{eff}}} + \rho \mathbf{g} + \mathbf{f} = \rho\ddot{\mathbf{u}},
 \f]
 * where \f$ \boldsymbol{\sigma_{\mathrm{eff}}} \f$ is the effective stress tensor,
 * \f$ \rho = (1 - \phi) \rho_s + \phi \rho_f \f$ is the average density of solids and fluids within the porous medium,
 * \f$ \mathbf{f} \f$ in \f$ \mathrm{N/m^3} \f$  is the external force acting on the body per unit volume (e.g. magnetism),
 * and \f$ \mathbf{u} = \mathbf{x} - \mathbf{x}_{\mathrm{initial}} \f$ is the displacement,
 * defined as the difference in material points \f$ \mathbf{x} \f$ and \f$ \mathbf{x}_{\mathrm{initial}} \f$
 * in the deformed and undeformed (initial) state, respectively. The model assumes quasi-static conditions,
 * that is, the above momentum balance equation is solved under the assumption that the acceleration term
 * \f$ \rho\ddot{\mathbf{u}} \approx 0\f$.
 *
 * Using the concept of the effective stress, the effective stress tensor \f$ \boldsymbol{\sigma_{\mathrm{eff}}} \f$ is
 * determined by the stress tensor \f$ \boldsymbol{\sigma} \f$ , the effective pore pressure \f$ p_{\mathrm{eff}} \f$ and the Biot's coefficient \f$ \alpha \f$ :
 \f[
 \boldsymbol{\sigma_{\mathrm{eff}}} = \boldsymbol{\sigma} - \alpha p_{\mathrm{eff}} \mathbf{I}
 \f]
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
#ifndef DUMUX_GEOMECHANICS_POROELASTIC_MODEL_HH
#define DUMUX_GEOMECHANICS_POROELASTIC_MODEL_HH

#include <dumux/poromechanics/poroelastic/model.hh>
#warning "This header is deprecated and will be removed after 3.10. Use dumux/poromechanics/poroelastic/model.hh."

 #endif
