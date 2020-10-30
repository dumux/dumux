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
 * \brief Calculates phase state for a single phase but two-component state.
 */
#ifndef DUMUX_PSEUDO1P2C_FLUID_STATE_HH
#define DUMUX_PSEUDO1P2C_FLUID_STATE_HH

#include <cassert>

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \brief Container for compositional variables in a 1p2c situation
 *
 *  This class holds variables for single-phase situations in a 2p2c context.
 *  It is used in case of a multiphysics approach. For the non-present phase,
 *  no information is stored but 0-values are returned to allow for general output
 *  methods.
 *  The "flash" calculation routines are in the sequential flash constrain solver, see
 *  CompositionalFlash .
 */
template <class ScalarType, class FluidSystem>
class PseudoOnePTwoCFluidState
{

public:
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    //! export the scalar type
    using Scalar = ScalarType;

    enum {
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx,

        comp0Idx = FluidSystem::comp0Idx,
        comp1Idx = FluidSystem::comp1Idx
    };

    /*! \name Acess functions */
    //@{
    /*!
     * \brief Returns the saturation \f$S_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * The saturation is defined as the pore space occupied by the fluid divided by the total pore space:
     *  \f[S_\alpha := \frac{\phi \mathcal{V}_\alpha}{\phi \mathcal{V}}\f]
     * This is set either to 1 or 0 depending on the phase presence.
     * \param phaseIdx the index of the phase
     */
    Scalar saturation(int phaseIdx) const
    { return phaseIdx == presentPhaseIdx_ ? 1.0 : 0.0; }

    //! \brief Returns the index of the phase that is present in that cell.
    int presentPhaseIdx() const
    { return presentPhaseIdx_; }

    /*!
     * \brief Return the partial pressure of a component in the gas phase.
     * \todo is this needed?
     *
     * For an ideal gas, this means \f$ R*T*c \f$.
     * Unit: \f$\mathrm{[Pa] = [N/m^2]}\f$
     *
     * \param compIdx the index of the component
     */
    Scalar partialPressure(int compIdx) const
    { return partialPressure(phase1Idx, compIdx); }

    /*!
     * \brief The partial pressure of a component in a phase \f$\mathrm{[Pa]}\f$
     * \todo is this needed?
     */
    Scalar partialPressure(int phaseIdx, int compIdx) const
    {
        assert(FluidSystem::isGas(phaseIdx));
        return pressure(phaseIdx)*moleFraction(phaseIdx, compIdx);
    }

    /*!
     * \brief The pressure \f$p_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa]}\f$
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief Set the density of a phase \f$\mathrm{[kg / m^3]}\f$
     */
    Scalar density(int phaseIdx) const
    { return phaseIdx == presentPhaseIdx_ ? density_ : 0.0; }

    /*!
     *  @copydoc CompositionalFluidState::molarDensity()
     */
    Scalar molarDensity(int phaseIdx) const
    { return phaseIdx == presentPhaseIdx_ ? molarDensity_ : 0.0; }

    /*!
     *  @copydoc CompositionalFluidState::massFraction()
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        if (phaseIdx != presentPhaseIdx_)
            return phaseIdx == compIdx ? 1.0 : 0.0;

        return compIdx == phase0Idx ? massFractionWater_ : 1.0 - massFractionWater_;
    }

    /*!
     * \brief Returns the molar fraction \f$x^\kappa_\alpha\f$ of the component \f$\kappa\f$ in fluid phase \f$\alpha\f$ in \f$\mathrm{[-]}\f$.
     *
     * This is either set to 1 or 0 depending on the phase presence for the
     * nonwetting phase in general.
     * It is set to the mole fraction of water or 1-moleFractionWater
     * if the considered component is the main component of the wetting phase.
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        if (phaseIdx != presentPhaseIdx_)
            return phaseIdx == compIdx ? 1.0 : 0.0;

        return compIdx == phase0Idx ? moleFractionWater_ : 1.0 - moleFractionWater_;
    }

    /*!
     * \brief The dynamic viscosity \f$\mu_\alpha\f$ of fluid phase \f$\alpha\f$ in \f$\mathrm{[Pa s]}\f$
     */
    Scalar viscosity(int phaseIdx) const
    {
        assert(phaseIdx == presentPhaseIdx_);
        return viscosity_;
    }

    /*!
     * \brief The average molar mass \f$\overline M_\alpha\f$ of phase \f$\alpha\f$ in \f$\mathrm{[kg/mol]}\f$
     *
     * The average molar mass is the mean mass of a mole of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \f[\mathrm{ \overline M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa}\f]
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return averageMolarMass_; }

    /*!
     * \brief The specific enthalpy \f$h_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     */
    Scalar enthalpy(int phaseIdx) const
    { return phaseIdx == presentPhaseIdx_ ? enthalpy_ : 0.0; }

    /*!
     * \brief The specific internal energy \f$u_\alpha\f$ of a fluid phase \f$\alpha\f$ in \f$\mathrm{[J/kg]}\f$
     *
     * The specific internal energy is defined by the relation:
     *
     * \f[u_\alpha = h_\alpha - \frac{p_\alpha}{\rho_\alpha}\f]
     */
    Scalar internalEnergy(int phaseIdx) const
    { return phaseIdx == presentPhaseIdx_ ?  enthalpy_ - pressure(phaseIdx)/density(phaseIdx) : 0.0; }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; }
    //@}

    /*!
     * \name Functions to set Data
     */
    //@{
    /*!
     * \brief Sets the viscosity of a phase \f$\mathrm{[Pa*s]}\f$.
     *
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setViscosity(int phaseIdx, Scalar value)
    {
        assert(phaseIdx == presentPhaseIdx_);
        viscosity_ = value;
    }

    /*!
     * \brief Sets the mass fraction of a component in a phase.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     * @param value Value to be stored
     */
    void setMassFraction(int phaseIdx, int compIdx, Scalar value)
    { massFractionWater_ = compIdx == comp0Idx ? value : 1.0 - value; }

    /*!
     * \brief Sets the molar fraction of a component in a fluid phase.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     * @param value Value to be stored
     */
    void setMoleFraction(int phaseIdx, int compIdx, Scalar value)
    { moleFractionWater_ = compIdx == comp0Idx ? value : 1.0 - value; }

    /*!
     * \brief Sets the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setDensity(int phaseIdx, Scalar value)
    {
        assert(phaseIdx == presentPhaseIdx_);
        density_ = value;
    }

    /*!
     * \brief Set the molar density of a phase \f$\mathrm{[mol / m^3]}\f$
     *
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setMolarDensity(int phaseIdx, Scalar value)
    {
        assert(phaseIdx == presentPhaseIdx_);
        molarDensity_ = value;
    }

    /*!
     * \brief Sets the phase Index that is present in this fluidState.
     * @param phaseIdx the index of the phase
     */
    void setPresentPhaseIdx(int phaseIdx)
    { presentPhaseIdx_ = phaseIdx; }

    /*!
     * \brief Sets the temperature
     *
     * @param value Value to be stored
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

    /*!
     * \brief Set the average molar mass of a fluid phase [kg/mol]
     *
     * The average molar mass is the mean mass of a mole of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \f[ \bar M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa \f]
     */
    void setAverageMolarMass(int phaseIdx, Scalar value)
    { averageMolarMass_ = value; }

    /*!
     * \brief Sets the phase pressure \f$\mathrm{[Pa]}\f$.
     */
    void setPressure(int phaseIdx, Scalar value)
    { pressure_[phaseIdx] = value; }

    /*!
     * \brief Sets phase enthalpy
     *
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setEnthalpy(int phaseIdx, Scalar value)
    {
        assert(phaseIdx == presentPhaseIdx_);
        enthalpy_ = value;
    }
    //@}

protected:
    Scalar pressure_[numPhases] = {};
    Scalar massConcentration_[numComponents] = {};
    Scalar averageMolarMass_ = 0.0;
    Scalar massFractionWater_ = 0.0;
    Scalar moleFractionWater_ = 0.0;
    Scalar density_ = 0.0;
    Scalar molarDensity_ = 0.0;
    Scalar viscosity_  = 0.0;
    Scalar enthalpy_ = 0.0;
    Scalar temperature_ = 0.0;
    int presentPhaseIdx_ = 0;
};

} // end namespace Dumux

#endif
