// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NIModel
 * \brief The implicit non-isothermal model.
 *
 * This model implements a generic energy balance for single and multi-phase
 * transport problems. Currently the non-isothermal model can be used on top of
 * the 1p2c, 2p, 2p2c and 3p3c models. Comparison to simple analytical solutions
 * for pure convective and conductive problems are found in the 1p2c test. Also refer
 * to this test for details on how to activate the non-isothermal model.
 *
 * For the energy balance, local thermal equilibrium is assumed. This
 * results in one energy conservation equation for the porous solid
 * matrix and the fluids:
 \f{align*}{
 \phi \frac{\partial \sum_\alpha \varrho_\alpha u_\alpha S_\alpha}{\partial t}
 & +
 \left( 1 - \phi \right) \frac{\partial (\varrho_s c_s T)}{\partial t}
 -
 \sum_\alpha \text{div}
 \left\{
 \varrho_\alpha h_\alpha
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left( \textbf{grad}\,p_\alpha - \varrho_\alpha \mathbf{g} \right)
 \right\} \\
    & - \text{div} \left(\lambda_{pm} \textbf{grad} \, T \right)
    - q^h = 0.
 \f}
 * <details>
 * <summary>Click for variables explanations:</summary>
 * \f$ \cdot \f$ \f$ \phi \f$ is the porosity of the porous-medium, <br>
 * \f$ \cdot \f$ \f$ S_\alpha \f$ represents the saturation of phase \f$ \alpha \f$, <br>
 * \f$ \cdot \f$ \f$ \rho_\alpha \f$ is the mass density of phase \f$ \alpha \f$, <br>
 * \f$ \cdot \f$ \f$ h_\alpha \f$ is the specific enthalpy of phase  \f$ \alpha \f$, <br>
 * \f$ \cdot \f$ \f$ u_\alpha \f$ is the specific internal energy of phase \f$ \alpha \f$, <br>
 * \f$ \cdot \f$ \f$ \lambda_{pm}\f$ is the heat conductivity in the porous medium, <br>
 * \f$ \cdot \f$ \f$ T \f$ is the Temperature, <br>
 * \f$ \cdot \f$ \f$ \rho_s \f$ is the mass density of phase the solid, <br>
 * \f$ \cdot \f$ \f$ c_s \f$ is the heat capacity of the solid, <br>
 * \f$ \cdot \f$ \f$ k_{r\alpha} \f$ is the relative permeability of phase \f$ \alpha \f$, <br>
 * \f$ \cdot \f$ \f$ \mu_\alpha \f$ is the dynamic viscosity of phase \f$ \alpha \f$, <br>
 * \f$ \cdot \f$ \f$ \mathbf{K} \f$ is the intrinsic permeability tensor, <br>
 * \f$ \cdot \f$ \f$ p_\alpha \f$ is the pressure of phase \f$ \alpha \f$, <br>
 * \f$ \cdot \f$ \f$ \mathbf{g} \f$ is the gravitational acceleration vector, <br>
 * \f$ \cdot \f$ \f$ q^h \f$ is a source or sink term.
 * </details>
 */

#ifndef DUMUX_NONISOTHERMAL_MODEL_HH
#define DUMUX_NONISOTHERMAL_MODEL_HH

#include <string>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>

namespace Dumux {

/*!
 * \ingroup NIModel
 * \brief Specifies a number properties of non-isothermal porous medium
 *        flow models based on the specifics of a given isothermal model.
 *
 * \tparam IsothermalTraits Model traits of the isothermal model
 */
template<class IsothermalT, class TDM = void>
struct PorousMediumFlowNIModelTraits : public IsothermalT
{
    //! Export the isothermal model traits
    using IsothermalTraits = IsothermalT;

    //! Export the thermal dispersion tensor type
    using ThermalDispersionModel = TDM;

    //! We solve for one more equation, i.e. the energy balance
    static constexpr int numEq() { return IsothermalTraits::numEq()+1; }
    //! only one energy equation is needed when assuming thermal equilibrium
    static constexpr int numEnergyEq() { return 1; }
    //! We additionally solve for the equation balance
    static constexpr bool enableEnergyBalance() { return true; }
    //! The indices related to the non-isothermal model
    using Indices = EnergyIndices< typename IsothermalTraits::Indices, numEq()>;
};

} // end namespace Dumux

#endif
