// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-     by Holger Class                                 *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Adaption of the BOX scheme to the non-isothermal three-phase three-component flow model.
 */
#ifndef DUMUX_NEW_3P3CNI_MODEL_HH
#define DUMUX_NEW_3P3CNI_MODEL_HH

#include <dumux/boxmodels/3p3c/3p3cmodel.hh>

namespace Dumux {
/*!
 * \ingroup ThreePThreeCNIModel
 * \brief Adaption of the BOX scheme to the non-isothermal three-phase three-component flow model.
 *
 * This model implements three-phase three-component flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$ each composed of up to three components
 * \f$\kappa \in \{ water, air, contaminant \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one transport equation for each component is obtained as
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_{\text{mol}, \alpha} x_\alpha^\kappa
 S_\alpha )}{\partial t}
 - \sum\limits_\alpha \nabla \cdot \left\{ \frac{k_{r\alpha}}{\mu_\alpha}
 \varrho_{\text{mol}, \alpha} x_\alpha^\kappa \mbox{\bf K}
 (\nabla p_\alpha - \varrho_{\text{mass}, \alpha} \mbox{\bf g}) \right\}
 \nonumber \\
 \nonumber \\
 && - \sum\limits_\alpha \nabla \cdot \left\{ D_{pm}^\kappa \varrho_{\text{mol},
 \alpha } \nabla x_\alpha^\kappa \right\}
 - q^\kappa = 0 \qquad \forall \kappa , \; \forall \alpha
 \f}
 *
 * Note that these balance equations above are molar.
 * In addition to that, a single balance of thermal energy is formulated
 * for the fluid-filled porous medium under the assumption of local thermal
 * equilibrium
 * \f{eqnarray*}
 && \phi \frac{\partial \left( \sum_\alpha \varrho_\alpha u_\alpha S_\alpha \right)}{\partial t}
 + \left( 1 - \phi \right) \frac{\partial (\varrho_s c_s T)}{\partial t}
 - \sum_\alpha \text{div} \left\{ \varrho_\alpha h_\alpha
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left( \text{grad}\,
 p_\alpha
 - \varrho_\alpha \mathbf{g} \right) \right\} \\
 &-& \text{div} \left( \lambda_{pm} \text{grad} \, T \right)
 - q^h = 0 \qquad \alpha \in \{w, n, g\}
 \f}
 *

 *
 * The equations are discretized using a fully-coupled vertex
 * centered finite volume (BOX) scheme as spatial scheme and
 * the implicit Euler method as temporal discretization.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations and
 * \f$x^w_\alpha + x^a_\alpha + x^c_\alpha = 1\f$ for the mole fractions.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are dependent on the locally present fluid phases
 * An adaptive primary variable switch is included. The phase state is stored for all nodes
 * of the system. The following cases can be distinguished:
 * <ul>
 *  <li> All three phases are present: Primary variables are two saturations (\f$S_w\f$ and \f$S_n\f$, a pressure, in this case \f$p_g\f$, and the temperature \f$T\f$. </li>
 *  <li> Only the water phase is present: Primary variables are now the mole fractions of air and contaminant in the water phase (\f$x_w^a\f$ and \f$x_w^c\f$), as well as temperature and the gas pressure, which is, of course, in a case where only the water phase is present, just the same as the water pressure. </li>
 *  <li> Gas and NAPL phases are present: Primary variables (\f$S_n\f$, \f$x_g^w\f$, \f$p_g\f$, \f$T\f$). </li>
 *  <li> Water and NAPL phases are present: Primary variables (\f$S_n\f$, \f$x_w^a\f$, \f$p_g\f$, \f$T\f$). </li>
 *  <li> Only gas phase is present: Primary variables (\f$x_g^w\f$, \f$x_g^c\f$, \f$p_g\f$, \f$T\f$). </li>
 *  <li> Water and gas phases are present: Primary variables (\f$S_w\f$, \f$x_w^g\f$, \f$p_g\f$, \f$T\f$). </li>
 * </ul>
 *
 */
template<class TypeTag>
class ThreePThreeCNIModel : public ThreePThreeCModel<TypeTag>
{
};

}

#include "3p3cnipropertydefaults.hh"

#endif
