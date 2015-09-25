// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief Adaption of the BOX scheme to the non-isothermal two-phase two-component flow model.
 */
#ifndef DUMUX_NEW_2P2CNI_MODEL_HH
#define DUMUX_NEW_2P2CNI_MODEL_HH

#include <dumux/boxmodels/2p2c/2p2cmodel.hh>

namespace Dumux {
/*!
 * \ingroup BoxModels
 * \defgroup TwoPTwoCNIModel Non-isothermal two-phase two-component box model
 */

/*!
 * \ingroup TwoPTwoCNIModel
 * \brief Adaption of the BOX scheme to the non-isothermal two-phase two-component flow model.
 *
 * This model implements a non-isothermal two-phase flow of two compressible and partly miscible fluids
 * \f$\alpha \in \{ w, n \}\f$. Thus each component \f$\kappa \{ w, a \}\f$ can be present in
 * each phase.
 * Using the standard multiphase Darcy approach a mass balance equation is
 * solved:
 * \f{eqnarray*}
    && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}{\partial t}
    - \sum_\alpha \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
    \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
    (\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g}) \right\}\\
    &-& \sum_\alpha \text{div} \left\{{\bf D_{\alpha, pm}^\kappa} \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
    - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
    \alpha \in \{w, n\}
 *     \f}
 * For the energy balance, local thermal equilibrium is assumed which results in one
 * energy conservation equation for the porous solid matrix and the fluids:
 * \f{eqnarray*}
    && \phi \frac{\partial \left( \sum_\alpha \varrho_\alpha u_\alpha S_\alpha \right)}{\partial t}
    + \left( 1 - \phi \right) \frac{\partial (\varrho_s c_s T)}{\partial t}
    - \sum_\alpha \text{div} \left\{ \varrho_\alpha h_\alpha
    \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left( \text{grad}\,
     p_\alpha
    - \varrho_\alpha \mathbf{g} \right) \right\} \\
    &-& \text{div} \left( \lambda_{pm} \text{grad} \, T \right)
    - q^h = 0 \qquad \alpha \in \{w, n\}
\f}
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * If both phases are present the primary variables are, like in the nonisothermal two-phase model, either \f$p_w\f$, \f$S_n\f$ and
 * temperature or \f$p_n\f$, \f$S_w\f$ and temperature. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPTwoIndices::pWsN</tt> or <tt>TwoPTwoCIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * In case that only one phase (nonwetting or wetting phase) is present the second primary
 * variable represents a mass fraction. The correct assignment of the second
 * primary variable is performed by a phase state dependent primary variable switch. The phase state is stored for all nodes of the system. The following cases can be distinguished:
 * <ul>
 *  <li>
 *    Both phases are present: The saturation is used (either\f$S_n\f$ or \f$S_w\f$, dependent on the chosen formulation).
 *  </li>
 *  <li>
 *    Only wetting phase is present: The mass fraction of air in the wetting phase \f$X^a_w\f$ is used.
 *  </li>
 *  <li>
 *    Only non-wetting phase is present: The mass fraction of water in the non-wetting phase, \f$X^w_n\f$, is used.
 *  </li>
 * </ul>
 */
template<class TypeTag>
class TwoPTwoCNIModel : public TwoPTwoCModel<TypeTag>
{
};

}

#include "2p2cnipropertydefaults.hh"

#endif
