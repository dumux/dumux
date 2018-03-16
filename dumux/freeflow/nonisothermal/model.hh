// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NavierStokesNIModel
 *
 * \brief A single-phase, non-isothermal Navier-Stokes model
 *
 * In addition to the momentum and mass/mole balance equations, this model also solves the <B> energy balance equation </B>:
 * \f[
 *    \frac{\partial (\varrho  v)}{\partial t}
 *    + \nabla \cdot \left( \varrho h {\boldsymbol{v}}
 *    - \lambda \textbf{grad}\, T \right) - q_T = 0
 * \f]
 *
 */

#ifndef DUMUX_STAGGERED_NI_MODEL_HH
#define DUMUX_STAGGERED_NI_MODEL_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesNIModel
 * \brief Specifies a number properties of non-isothermal free flow
 *        flow models based on the specifics of a given isothermal model.
 * \tparam IsothermalTraits Model traits of the isothermal model
 */
template<class IsothermalTraits>
struct NavierStokesNIModelTraits : public IsothermalTraits
{
    //! We solve for one more equation, i.e. the energy balance
    static constexpr int numEq() { return IsothermalTraits::numEq()+1; }
    //! We additionally solve for the equation balance
    static constexpr bool enableEnergyBalance() { return true; }
};

} // end  namespace Dumux

#endif
