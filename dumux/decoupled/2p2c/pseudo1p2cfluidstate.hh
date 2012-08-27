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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Calculates phase state for a single phase but two-component state.
 */
#ifndef DUMUX_PSEUDO1P2C_FLUID_STATE_HH
#define DUMUX_PSEUDO1P2C_FLUID_STATE_HH

#include <dumux/decoupled/2p2c/2p2cproperties.hh>

namespace Dumux
{
/*!
 * \ingroup multiphysics
 * \brief Container for compositional variables in a 1p2c situation
 *
 *  This class represents a pseudo flash calculation when only single-phase situations
 *  occur. This is used in case of a multiphysics approach.
 *  \tparam TypeTag The property Type Tag
 */
template <class TypeTag>
class PseudoOnePTwoCFluidState
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

public:
    enum {     numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
            numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = Indices::wPhaseIdx,
        nCompIdx = Indices::nPhaseIdx
    };

public:
    /*!
     * \name flash calculation routines
     * Routine to determine the phase composition after the transport step.
     */
    //@{
    //! The simplest possible update routine for 1p2c "flash" calculations
    /*!
     * Routine goes as follows:
     * - Check if we are in single phase condition
     * - Assign total concentration to the present phase
     *
     * \param Z1 Feed mass fraction \f$\mathrm{[-]}\f$
     * \param phasePressure Vector holding the pressure \f$\mathrm{[Pa]}\f$
     * \param presentPhaseIdx Subdomain Index = Indication which phase is present
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    /*! \name Acess functions */
    //@{
    /*! \brief Returns the saturation of a phase.
     *  \param phaseIdx Index of the phase
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == presentPhaseIdx_)
            return 1.;
        else
            return 0.;
    };

    int presentPhaseIdx() const
    {
        return presentPhaseIdx_;
    }

    /*! \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *  \param phaseIdx Index of the phase
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    Scalar density(int phaseIdx) const
    {
        if(phaseIdx == presentPhaseIdx_)
            return density_;
        else
            return 0.; // only required for output if phase is not present anyhow
    }

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     *  \param phaseIdx Index of the phase
     *  \param compIdx Index of the component
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        if(phaseIdx != presentPhaseIdx_)
            return 0.;


        if (compIdx == wPhaseIdx)
            return massFractionWater_[phaseIdx];
        else
            return 1.-massFractionWater_[phaseIdx];

    }

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     *  \param phaseIdx Index of the phase
     *  \param compIdx Index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        if(phaseIdx != presentPhaseIdx_)
            return 0.;
        if (compIdx == wPhaseIdx)
            return moleFractionWater_[phaseIdx];
        else
            return 1.-moleFractionWater_[phaseIdx];
    }

    /*!
     * \brief Returns the viscosity of a phase TODO: \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar viscosity(int phaseIdx) const
    {
        assert(phaseIdx == presentPhaseIdx_);
        return viscosity_;
    }

    Scalar averageMolarMass(int phaseIdx) const
    {
        return aveMoMass_;
    }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; };
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
    {
        if (compIdx == wCompIdx)
            massFractionWater_[phaseIdx] = value;
        else
            massFractionWater_[phaseIdx] = 1- value;
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
        if (compIdx == wCompIdx)
            moleFractionWater_[phaseIdx] = value;
        else
            moleFractionWater_[phaseIdx] = 1-value;
    }
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
     * \brief Sets the phase Index that is present in this fluidState.
     * @param phaseIdx the index of the phase
     */
    void setPresentPhaseIdx(int phaseIdx)
    {
        presentPhaseIdx_ = phaseIdx;
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
    //TODO: doc me
    void setAverageMolarMass(int phaseIdx, Scalar value)
    {
        aveMoMass_ = value;
    }
    /*!
     * \brief Sets the phase pressure \f$\mathrm{[Pa]}\f$.
     */
    void setPressure(int phaseIdx, Scalar value)
    {
        pressure_[phaseIdx] = value;
    };
    //@}

public:
    Scalar aveMoMass_;
    Scalar massConcentration_[numComponents];
    Scalar massFractionWater_[numPhases];
    Scalar moleFractionWater_[numPhases];
    Scalar pressure_[numPhases];
    Scalar density_;
    Scalar viscosity_;
    Scalar temperature_;
    int presentPhaseIdx_;
};

} // end namepace

#endif
