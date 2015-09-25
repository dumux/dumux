// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
 * \brief Determines the phase state of the immiscible sequential 2p model.
 */
#ifndef DUMUX_2P_FLUID_STATE_HH
#define DUMUX_2P_FLUID_STATE_HH

#include <dumux/material/fluidstate.hh>
#include <dumux/decoupled/2p/2pproperties.hh>

namespace Dumux
{
/*!
 * \ingroup IMPES
 */
/*!
 * \brief Calcultes the phase state from the primary variables in the sequential
 *        2p model.
 */
template <class TypeTag>
class TwoPFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                           TwoPFluidState<TypeTag> >
{
    typedef TwoPFluidState<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = Indices::wPhaseIdx,
        nCompIdx = Indices::nPhaseIdx,
    };

public:
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };

public:
    void update(Scalar Sw, Scalar pressW, Scalar pressN, Scalar temperature)
    {
        Sw_ = Sw;
        phasePressure_[wPhaseIdx] = pressW;
        phasePressure_[nPhaseIdx] = pressN;
        temperature_=temperature;
        density_[wPhaseIdx] = FluidSystem::phaseDensity(wPhaseIdx,
                                                        temperature,
                                                        pressW,
                                                        *this);
        density_[nPhaseIdx] = FluidSystem::phaseDensity(nPhaseIdx,
                                                        temperature,
                                                        pressN,
                                                        *this);
    }

    /*!
     * \brief Returns the saturation of a phase.
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
            return Sw_;
        else
            return Scalar(1.0) - Sw_;
    };

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    {
        if (compIdx == phaseIdx)
            return 1.0;
        return 0;
    }

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    {
        return massFrac(phaseIdx, compIdx);
    }

    /*!
     * \brief Returns the total concentration of a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * This is equivalent to the sum of all component concentrations.
     */
    Scalar phaseConcentration(int phaseIdx) const
    {
        return density_[phaseIdx]/FluidSystem::molarMass(phaseIdx);
    };

    /*!
     * \brief Returns the concentration of a component in a phase \f$\mathrm{[mol/m^3]}\f$.
     */
    Scalar concentration(int phaseIdx, int compIdx) const
    {
        if (phaseIdx == compIdx)
            return phaseConcentration(phaseIdx);
        return 0;
    };

    /*!
     * \brief Returns the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Returns mean molar mass of a phase \f$\mathrm{[kg/mol]}\f$.
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return FluidSystem::molarMass(phaseIdx); };

    /*!
     * \brief Returns the partial pressure of a component in the gas phase \f$\mathrm{[Pa]}\f$.
     */
    Scalar partialPressure(int compIdx) const
    {
        if (compIdx == wPhaseIdx)
            return 0;
        return phasePressure_[nPhaseIdx];
    }

    /*!
     * \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     */
    Scalar phasePressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure \f$\mathrm{[Pa]}\f$
     */
    Scalar capillaryPressure() const
    { return phasePressure_[nPhaseIdx] - phasePressure_[wPhaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

public:
    Scalar density_[numPhases];
    Scalar phasePressure_[numPhases];
    Scalar temperature_;
    Scalar Sw_;
};

} // end namepace

#endif
