// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \brief Adaption of the BOX scheme to the non-isothermal twophase flow model.
 */

#ifndef DUMUX_2PNI_MODEL_HH
#define DUMUX_2PNI_MODEL_HH

#include <dumux/boxmodels/2p/2pmodel.hh>

namespace Dumux {

/*!
 * \ingroup BoxModels
 * \defgroup TwoPNIBoxModel Non-isothermal two-phase box model
 */

/*!
 * \ingroup TwoPNIBoxModel
 * \brief Adaption of the BOX scheme to the non-isothermal twophase flow model.
 *
 * This model implements a non-isothermal two-phase flow of two completely immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$.
 * Using the standard multiphase Darcy approach, the mass conservation equations for both phases can
 * be described as follows:
 * \f{eqnarray*}
 && \phi \frac{\partial (\varrho_\alpha S_\alpha )}{\partial t}
 - \text{div} \left\{ \varrho_\alpha
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g}) \right\}
 - q_\alpha^\kappa = 0 \qquad \alpha \in \{w, n\}
 *     \f}
 * For the energy balance, local thermal equilibrium is assumed which results in one
 * energy conservation equation for the porous solid matrix and the fluids:
 * \f{eqnarray*}
 && \phi \frac{\partial \left( \sum_\alpha \varrho_\alpha u_\alpha S_\alpha \right)}{\partial t}
 + \left( 1 - \phi \right) \frac{\partial (\varrho_s c_s T)}{\partial t}
 - \sum_\alpha \text{div} \left\{ \varrho_\alpha h_\alpha
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K} \left( \text{grad} \, p_\alpha
 - \varrho_\alpha \mbox{\bf g} \right) \right\} \\
    &-& \text{div} \left( \lambda_{pm} \text{grad} \, T \right)
 - q^h = 0, \qquad \alpha \in \{w, n\}.
 \f}
 *
 * The equations are discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and the implicit Euler method
 * as time discretization.
 *
 * Currently the model supports choosing either \f$p_w\f$, \f$S_n\f$ and \f$T\f$ or \f$p_n\f$,
 * \f$S_w\f$ and \f$T\f$ as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either <tt>TwoPNIIndices::pWsN</tt> or <tt>TwoPIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$, \f$S_n\f$ and \f$T\f$.
 */
template<class TypeTag>
class TwoPNIModel: public TwoPModel<TypeTag>
{
};

}

#include "2pnipropertydefaults.hh"

#endif
