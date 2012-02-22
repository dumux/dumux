// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \brief Adaption of the BOX scheme to the non-isothermal twophase flow model.
 */

#ifndef DUMUX_2PNI_MODEL_HH
#define DUMUX_2PNI_MODEL_HH

#include <dumux/boxmodels/2p/2pmodel.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNIBoxModel
 * \brief A two-phase, non-isothermal flow model using the box scheme.
 *
 * This model implements a non-isothermal two-phase flow for two
 * immiscible fluids \f$\alpha \in \{ w, n \}\f$. Using the standard
 * multiphase Darcy approach, the mass conservation equations for both
 * phases can be described as follows:
 * \f[
 \phi \frac{\partial \phi \varrho_\alpha S_\alpha}{\partial t}
 - 
 \text{div} 
 \left\{ 
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathrm{K}
 \left( \textrm{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 \right\}
 -
 q_\alpha = 0 \qquad \alpha \in \{w, n\}
 \f]
 *
 * For the energy balance, local thermal equilibrium is assumed. This
 * results in one energy conservation equation for the porous solid
 * matrix and the fluids: 
 
 \f{eqnarray*}
 \frac{\partial \phi \sum_alpha \varrho_\alpha u_\alpha S_\alpha}{\partial t}
 & + 
 \left( 1 - \phi \right) \frac{\partial (\varrho_s c_s T)}{\partial t}
 - 
 \sum_\alpha \text{div} 
 \left\{
 \varrho_\alpha h_\alpha
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} 
 \left( \textbf{grad}\,p_\alpha - \varrho_\alpha \mbox{\bf g} \right)
 \right\} \                                                     \
    & - \text{div} \left(\lambda_{pm} \textbf{grad} \, T \right)
    - q^h = 0, \qquad \alpha \in \{w, n\} \;,
 \f}
 * where \f$h_\alpha\f$ is the specific enthalpy of a fluid phase
 * \f$\alpha\f$ and \f$u_\alpha = h_\alpha -
 * p_\alpha/\varrho_\alpha\f$ is the specific internal energy of the
 * phase.
 *
 * The equations are discretized using a fully-coupled vertex centered
 * finite volume (box) scheme as spatial and the implicit Euler method
 * as time discretization.
 *
 * Currently the model supports choosing either \f$p_w\f$, \f$S_n\f$
 * and \f$T\f$ or \f$p_n\f$, \f$S_w\f$ and \f$T\f$ as primary
 * variables. The formulation which ought to be used can be specified
 * by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPNIIndices::pWsN</tt> or <tt>TwoPIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$, \f$S_n\f$ and \f$T\f$.
 */
template<class TypeTag>
class TwoPNIModel: public TwoPModel<TypeTag>
{};

}

#include "2pnipropertydefaults.hh"

#endif
