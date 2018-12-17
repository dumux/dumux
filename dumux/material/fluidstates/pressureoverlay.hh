// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup FluidStates
 * \brief This is a fluid state which allows to set the fluid
 *        pressures and takes all other quantities from an other
 *        fluid state.
 */
#ifndef DUMUX_PRESSURE_OVERLAY_FLUID_STATE_HH
#define DUMUX_PRESSURE_OVERLAY_FLUID_STATE_HH

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \brief This is a fluid state which allows to set the fluid
 *        pressures and takes all other quantities from an other
 *        fluid state.
 */
template <class FluidState>
class PressureOverlayFluidState
{
public:
    static constexpr int numPhases = FluidState::numPhases;
    static constexpr int numComponents = FluidState::numComponents;

    //! export the scalar type
    using Scalar = typename FluidState::Scalar;

    /*!
     * \brief Constructor
     *
     * \param fs Fluidstate
     * The overlay fluid state copies the pressures from the argument,
     * so it initially behaves exactly like the underlying fluid
     * state.
     */
    PressureOverlayFluidState(const FluidState &fs)
    : fs_(&fs)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            pressure_[phaseIdx] = fs.pressure(phaseIdx);
    }

    // copy & move constructor / assignment operators
    PressureOverlayFluidState(const PressureOverlayFluidState &fs) = default;
    PressureOverlayFluidState(PressureOverlayFluidState &&fs) = default;
    PressureOverlayFluidState& operator=(const PressureOverlayFluidState &fs) = default;
    PressureOverlayFluidState& operator=(PressureOverlayFluidState &&fs) = default;

    /*****************************************************
     * Generic access to fluid properties (No assumptions
     * on thermodynamic equilibrium required)
     *****************************************************/
    /*!
     * \brief Returns the saturation \f$S_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The saturation is defined as the pore space occupied by the fluid divided by the total pore space:
     *  \f[S_\alpha := \frac{\phi \mathcal{V}_\alpha}{\phi \mathcal{V}}\f]
     *
     * \param phaseIdx the index of the phase
     */
    Scalar saturation(int phaseIdx) const
    { return fs_->saturation(phaseIdx); }

    /*!
     * \brief Returns the molar fraction \f$x^\kappa_\alpha\f$ of the component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The molar fraction \f$x^\kappa_\alpha\f$ is defined as the ratio of the number of molecules
     * of component \f$\kappa\f$ and the total number of molecules of the phase \f$\alpha\f$.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    { return fs_->moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mass fraction \f$X^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The mass fraction \f$X^\kappa_\alpha\f$ is defined as the weight of all molecules of a
     * component divided by the total mass of the fluid phase. It is related with the component's mole fraction by means of the relation
     *
     * \f[X^\kappa_\alpha = x^\kappa_\alpha \frac{M^\kappa}{\overline M_\alpha}\;,\f]
     *
     * where \f$M^\kappa\f$ is the molar mass of component \f$\kappa\f$ and
     * \f$\overline M_\alpha\f$ is the mean molar mass of a molecule of phase
     * \f$\alpha\f$.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    { return fs_->massFraction(phaseIdx, compIdx); }

    /*!
     * \brief The average molar mass \f$\overline M_\alpha\f$ of phase \f$\alpha\f$ in \f$\mathrm{[kg/mol]}\f$
     *
     * The average molar mass is the mean mass of a mole of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \f[\mathrm{ \overline M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa}\f]
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return fs_->averageMolarMass(phaseIdx); }

    /*!
     * \brief The molar concentration \f$c^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * This quantity is usually called "molar concentration" or just
     * "concentration", but there are many other (though less common)
     * measures for concentration.
     *
     * http://en.wikipedia.org/wiki/Concentration
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return fs_->molarity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity \f$f^\kappa_\alpha\f$ of component \f$\kappa\f$
     *  in fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa]}\f$
     *
     *  The fugacity is defined as:
     *  \f$f_\alpha^\kappa := \Phi^\kappa_\alpha x^\kappa_\alpha p_\alpha \;,\f$
     *  where \f$\Phi^\kappa_\alpha\f$ is the fugacity coefficient \cite reid1987 .
     *  The physical meaning of fugacity becomes clear from the equation:
     *       \f[f_\alpha^\kappa = p_\alpha \exp\left\{\frac{\zeta^\kappa_\alpha}{R T_\alpha} \right\} \;,\f]
     *  where \f$\zeta^\kappa_\alpha\f$ represents the \f$\kappa\f$'s chemical
     *  potential in phase \f$\alpha\f$, \f$R\f$ stands for the ideal gas constant,
     *  and \f$T_\alpha\f$ for the absolute temperature of phase \f$\alpha\f$. Assuming thermal equilibrium,
     *  there is a one-to-one mapping between a component's chemical potential
     *  \f$\zeta^\kappa_\alpha\f$ and its fugacity \f$f^\kappa_\alpha\f$. In this
     *  case chemical equilibrium can thus be expressed by:
     *     \f[f^\kappa := f^\kappa_\alpha = f^\kappa_\beta\quad\forall \alpha, \beta\f]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return fs_->fugacity(phaseIdx, compIdx); }

    /*!
     * \brief The fugacity coefficient \f$\Phi^\kappa_\alpha\f$ of component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fs_->fugacityCoefficient(phaseIdx, compIdx); }

    /*!
     * \brief The molar volume \f$v_{mol,\alpha}\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[m^3/mol]}\f$
     *
     * This quantity is the inverse of the molar density.
     */
    Scalar molarVolume(int phaseIdx) const
    { return fs_->molarVolume(phaseIdx); }

    /*!
     * \brief The mass density \f$\rho_\alpha\f$ of the fluid phase
     *  \f$\alpha\f$ in \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar density(int phaseIdx) const
    { return fs_->density(phaseIdx); }

    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    Scalar molarDensity(int phaseIdx) const
    { return fs_->molarDensity(phaseIdx); }

    /*!
     * \brief The absolute temperature\f$T_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[K]}\f$
     */
    Scalar temperature(int phaseIdx) const
    { return fs_->temperature(phaseIdx); }

    /*!
     * \brief The pressure \f$p_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa]}\f$
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief The specific enthalpy \f$h_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     */
    Scalar enthalpy(int phaseIdx) const
    { return fs_->enthalpy(phaseIdx); }

    /*!
     * \brief The specific internal energy \f$u_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     *
     * The specific internal energy is defined by the relation:
     *
     * \f[u_\alpha = h_\alpha - \frac{p_\alpha}{\rho_\alpha}\f]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return fs_->internalEnergy(phaseIdx); }

    /*!
     * \brief The dynamic viscosity \f$\mu_\alpha\f$ of fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa s]}\f$
     */
    Scalar viscosity(int phaseIdx) const
    { return fs_->viscosity(phaseIdx); }


    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the pressure \f$\mathrm{[Pa]}\f$ of a fluid phase
     */
    void setPressure(int phaseIdx, Scalar value)
    { pressure_[phaseIdx] = value; }

protected:
    const FluidState *fs_;
    Scalar pressure_[numPhases] = {};
};

} // end namespace Dumux

#endif
