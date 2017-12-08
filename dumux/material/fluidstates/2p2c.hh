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
 *
 * \brief Calculates the 2p2c phase state for compositional models.
 */
#ifndef DUMUX_2P2C_FLUID_STATE_HH
#define DUMUX_2P2C_FLUID_STATE_HH

#include <dumux/common/valgrind.hh>

namespace Dumux
{
/*!
 * \ingroup FluidStates
 * \brief Calculates the phase state from the primary variables in the
 *        sequential 2p2c model.
 *
 *        This boils down to so-called "flash calculation", in this case isothermal and isobaric.
 *
 *  \tparam TypeTag The property Type Tag
 */
template <class Scalar, class FluidSystem>
class TwoPTwoCFluidState
{
public:
    enum {
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
    };

    enum {  numPhases = FluidSystem::numPhases,
            numComponents = FluidSystem::numComponents};
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    /*!
     * \name acess functions
     */
    //@{
    /*! @copydoc CompositionalFluidState::saturation()
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
        {
            return sw_;
        }
        return 1.0 - sw_;
    }

    /*! @copydoc CompositionalFluidState::massFraction()
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        return massFraction_[phaseIdx][compIdx];
    }

    /*! @copydoc CompositionalFluidState::moleFraction()
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        return moleFraction_[phaseIdx][compIdx];
    }

    /*! @copydoc CompositionalFluidState::density()
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*! @copydoc CompositionalFluidState::viscosity()
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*! @copydoc CompositionalFluidState::partialPressure()
     */
    Scalar partialPressure(int compIdx) const
    {
        return partialPressure(nPhaseIdx, compIdx);
    }

    /*! @copydoc CompositionalFluidState::partialPressure()
     */
    Scalar partialPressure(int phaseIdx, int compIdx) const
    {
        assert(FluidSystem::isGas(phaseIdx));
        return pressure(phaseIdx)*moleFraction(phaseIdx, compIdx);
    }

    /*! @copydoc CompositionalFluidState::pressure()
     */
    Scalar pressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure \f$\mathrm{[Pa]}\f$
     */
    Scalar capillaryPressure() const
    { return phasePressure_[nPhaseIdx] - phasePressure_[wPhaseIdx]; }

    /*! @copydoc CompositionalFluidState::temperature()
     */
    Scalar temperature(int phaseIdx = 0) const
    { return temperature_; }

    /*! @copydoc CompositionalFluidState::averageMolarMass()
     */
    Scalar averageMolarMass(int phaseIdx) const
    {
        Scalar averageMolarMass = 0;

        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            averageMolarMass += moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);
        }
        return averageMolarMass;
    }

    /*!
     * \brief Returns the phase mass fraction. phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar phaseMassFraction(int phaseIdx)
    {
      using std::isnan;
        if (isnan(nu_[phaseIdx]))  //in contrast to the standard update() method, satflash() does not calculate nu.
        {
            nu_[wPhaseIdx] = sw_ * density_[wPhaseIdx] / (sw_ * density_[wPhaseIdx] + (1. - sw_) * density_[nPhaseIdx]);
            nu_[nPhaseIdx] = 1. - nu_[wPhaseIdx];
            return nu_[phaseIdx];
        }
        return nu_[phaseIdx];
    }
    /*!
     * \brief Returns the phase mass fraction \f$ \nu \f$:
     *  phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar& nu(int phaseIdx) const
    {
        return phaseMassFraction(phaseIdx);
    }

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
    { viscosity_[phaseIdx] = value; }


    /*!
     * \brief Sets the mass fraction of a component in a phase.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     * @param value Value to be stored
     */
    void setMassFraction(int phaseIdx, int compIdx, Scalar value)
    {
        massFraction_[phaseIdx][compIdx] = value;
    }

    /*!
     * \brief Sets the molar fraction of a component in a fluid phase.
     *
     * \param phaseIdx the index of the phase
     * \param compIdx the index of the component
     * @param value Value to be stored
     */
    void setMoleFraction(int phaseIdx, int compIdx, Scalar value)
    {
        moleFraction_[phaseIdx][compIdx] = value;
    }
    /*!
     * \brief Sets the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }
    /*!
     * \brief Sets the saturation of a phase.
     * Internally, only the wetting saturation is stored.
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setSaturation(int phaseIdx, Scalar value)
    {
        if (phaseIdx == wPhaseIdx)
            sw_ = value;
        else
            sw_ = 1.-value;
    }

    /*!
     * \brief Sets the phase mass fraction. phase mass per total mass \f$\mathrm{[kg/kg]}\f$.
     *
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setNu(int phaseIdx, Scalar value)
    {
            nu_[phaseIdx] = value;
    }
    /*!
     * \brief Sets the temperature
     *
     * @param value Value to be stored
     */
    void setTemperature(Scalar value)
    {
        temperature_ = value;
    }
    /*!
     * \brief Sets phase pressure
     *
     * \param phaseIdx the index of the phase
     * @param value Value to be stored
     */
    void setPressure(int phaseIdx, Scalar value)
    {
        phasePressure_[phaseIdx] = value;
    }

    //@}
    TwoPTwoCFluidState()
    { Valgrind::SetUndefined(*this); }

protected:
//    Scalar massConcentration_[numComponents];
    Scalar phasePressure_[numPhases];
    Scalar temperature_;

    Scalar sw_;
    PhaseVector nu_; //phase mass fraction
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar massFraction_[numPhases][numComponents];
    Scalar moleFraction_[numPhases][numComponents];
};

} // end namespace Dumux

#endif
