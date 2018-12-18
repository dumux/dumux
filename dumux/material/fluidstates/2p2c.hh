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
 * \brief Calculates the 2p2c phase state for compositional models.
 */
#ifndef DUMUX_2P2C_FLUID_STATE_HH
#define DUMUX_2P2C_FLUID_STATE_HH

#include <cmath>
#include <dune/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \brief Calculates the phase state from the primary variables in the
 *        sequential 2p2c model.
 * This boils down to so-called "flash calculation", in this case isothermal and isobaric.
 */
template <class ScalarType, class FluidSystem>
class DUNE_DEPRECATED_MSG("Use CompositionalFluidState instead!") TwoPTwoCFluidState
{
public:
    enum {
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx,
    };

public:
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    //! export the scalar type
    using Scalar = ScalarType;

    // comply with new style 2p2c models
    int wettingPhase() const
    { return phase0Idx; }

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
    { return phaseIdx == phase0Idx ? sw_ : 1.0 - sw_; }

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
    { return moleFraction_[phaseIdx][compIdx]; }

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
    { return massFraction_[phaseIdx][compIdx]; }

    /*!
     * \brief The mass density \f$\rho_\alpha\f$ of the fluid phase
     *  \f$\alpha\f$ in \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*! @copydoc CompositionalFluidState::molarDensity()
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; }

    /*! @copydoc CompositionalFluidState::viscosity()
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*!
     * \brief The partial pressure of a component in the n-phase \f$\mathrm{[Pa]}\f$
     * \todo is this necessary?
     */
    Scalar partialPressure(int compIdx) const
    { return partialPressure(phase1Idx, compIdx); }

    /*!
     * \brief The partial pressure of a component in a phase \f$\mathrm{[Pa]}\f$
     * \todo is this necessary?
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
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure \f$\mathrm{[Pa]}\f$
     */
    Scalar capillaryPressure() const
    { return phasePressure_[phase1Idx] - phasePressure_[phase0Idx]; }

    /*!
     * \brief The temperature within the domain \f$\mathrm{[K]}\f$
     */
    Scalar temperature(int phaseIdx = 0) const
    { return temperature_; }

    /*!
     * \brief The average molar mass \f$\overline M_\alpha\f$ of phase \f$\alpha\f$ in \f$\mathrm{[kg/mol]}\f$
     *
     * The average molar mass is the mean mass of a mole of the
     * fluid at current composition. It is defined as the sum of the
     * component's molar masses weighted by the current mole fraction:
     * \f[\mathrm{ \overline M_\alpha = \sum_\kappa M^\kappa x_\alpha^\kappa}\f]
     */
    Scalar averageMolarMass(int phaseIdx) const
    {
        Scalar averageMolarMass = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            averageMolarMass += moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);

        return averageMolarMass;
    }

    /*!
     * \brief Returns the phase mass fraction. phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     * \param phaseIdx the index of the phase
     */
    Scalar phaseMassFraction(int phaseIdx)
    {
        using std::isnan;
        if (isnan(nu_[phaseIdx]))  //in contrast to the standard update() method, satflash() does not calculate nu.
        {
            nu_[phase0Idx] = sw_ * density_[phase0Idx] / (sw_ * density_[phase0Idx] + (1. - sw_) * density_[phase1Idx]);
            nu_[phase1Idx] = 1. - nu_[phase0Idx];
            return nu_[phaseIdx];
        }
        return nu_[phaseIdx];
    }

    /*!
     * \brief Returns the phase mass fraction \f$ \nu \f$:
     *  phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     * \param phaseIdx the index of the phase
     */
    Scalar nu(int phaseIdx) const
    {
        return phaseMassFraction(phaseIdx);
    }

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Sets the viscosity of a phase \f$\mathrm{[Pa*s]}\f$.
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setViscosity(int phaseIdx, Scalar value)
    { viscosity_[phaseIdx] = value; }


    /*!
     * \brief Sets the mass fraction of a component in a phase.
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     * @param value Value to be stored
     */
    void setMassFraction(int phaseIdx, int compIdx, Scalar value)
    { massFraction_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Sets the molar fraction of a component in a fluid phase.
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     * @param value Value to be stored
     */
    void setMoleFraction(int phaseIdx, int compIdx, Scalar value)
    { moleFraction_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Sets the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Set the molar density of a phase \f$\mathrm{[mol / m^3]}\f$
     */
    void setMolarDensity(int phaseIdx, Scalar value)
    { molarDensity_[phaseIdx] = value; }

    /*!
     * \brief Sets the saturation of a phase.
     * Internally, only the wetting saturation is stored.
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setSaturation(int phaseIdx, Scalar value)
    { sw_ = phaseIdx == phase0Idx ? value : 1.0 - value; }

    /*!
     * \brief Sets the phase mass fraction. phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setNu(int phaseIdx, Scalar value)
    { nu_[phaseIdx] = value; }

    /*!
     * \brief Sets the temperature
     * @param value Value to be stored
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

    /*!
     * \brief Sets phase pressure
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setPressure(int phaseIdx, Scalar value)
    { phasePressure_[phaseIdx] = value; }

protected:
    //! zero-initialize all data members with braces syntax
    Scalar temperature_ = 0.0;
    Scalar sw_ = 0.0;
    Scalar phasePressure_[numPhases] = {};
    Scalar nu_[numPhases] = {};
    Scalar density_[numPhases] = {};
    Scalar molarDensity_[numPhases] = {};
    Scalar viscosity_[numPhases] = {};
    Scalar massFraction_[numPhases][numComponents] = {};
    Scalar moleFraction_[numPhases][numComponents] = {};
};

} // end namespace Dumux

#endif
