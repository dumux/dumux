// $Id: 1p2cfluidstate.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                               *
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
 * \brief Calcultes the phase state from the primary variables in the
 *        1p2c model.
 */
#ifndef DUMUX_1P2C_PHASE_STATE_HH
#define DUMUX_1P2C_PHASE_STATE_HH

#include "1p2cproperties.hh"

#include <dumux/material/fluidstate.hh>

namespace Dumux
{
/*!
 * \brief Calcultes the phase state from the primary variables in the
 *        1p2c model.
 */
template <class TypeTag>
class OnePTwoCFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                             OnePTwoCFluidState<TypeTag> >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum {
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,

        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx,

        phaseIndex = GET_PROP_VALUE(TypeTag, PTAG(PhaseIndex)),
        comp1Index = GET_PROP_VALUE(TypeTag, PTAG(Comp1Index)),
        comp2Index = GET_PROP_VALUE(TypeTag, PTAG(Comp2Index)),
    };

public:
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };

    /*!
     * \brief Update the phase state from the primary variables.
     *
     * \param primaryVars The primary variables
     * \param temperature The temperature in Kelvin
     */
    void update(const PrimaryVariables &primaryVars,
                Scalar temperature)
    {
        Valgrind::CheckDefined(primaryVars);

        temperature_ = temperature;

        phasePressure_ = primaryVars[pressureIdx];
        x1_ = primaryVars[x1Idx];
        meanMolarMass_ =
            (1 - x1_)*FluidSystem::molarMass(comp1Index) +
            (x1_    )*FluidSystem::molarMass(comp2Index);

        density_ = FluidSystem::phaseDensity(phaseIndex, temperature_, phasePressure_, *this);
        molarDensity_ = density_ / meanMolarMass_;

        Valgrind::CheckDefined(x1_);
        Valgrind::CheckDefined(phasePressure_);
        Valgrind::CheckDefined(density_);
        Valgrind::CheckDefined(meanMolarMass_);
        Valgrind::CheckDefined(temperature_);
        Valgrind::CheckDefined(*this);
    }

    /*!
     * \brief Returns the saturation of a phase.
     *
     * \param phaseIdx The index of the considered phase
     */
    Scalar saturation(int phaseIdx) const
    { return (phaseIndex == phaseIdx)?1.0:0.0; };

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     *
     * \param phaseIdx The index of the considered phase
     * \param compIdx The index of the considered component
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    {
        // we are a single phase model!
        if (phaseIdx != phaseIndex) return 0.0;

        if (compIdx==comp1Index)
            return 1-x1_;
        else if (compIdx==comp2Index)
            return x1_;
        return 0.0;
    }

    /*!
     * \brief Returns the total concentration of a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * This is equivalent to the sum of all component concentrations.
     * \param phaseIdx The index of the considered phase
     */
    Scalar phaseConcentration(int phaseIdx) const
    {
        if (phaseIdx != phaseIndex)
            return 0;
        return density_/meanMolarMass_;
    };

    /*!
     * \brief Returns the concentration of a component in a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param phaseIdx The index of the considered phase
     * \param compIdx The index of the considered component
     */
    Scalar concentration(int phaseIdx, int compIdx) const
    { return phaseConcentration(phaseIdx)*moleFrac(phaseIdx, compIdx); };


    /*!
     * \brief Returns the mass fraction of a component in a phase.
     *
     * \param phaseIdx The index of the considered phase
     * \param compIdx The index of the considered component
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    {
        if (phaseIdx != phaseIndex)
            return 0;
        return
            moleFrac(phaseIdx, compIdx)*
            FluidSystem::molarMass(compIdx)
            /
            meanMolarMass_;
    }

    /*!
     * \brief Returns the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx The index of the considered phase
     *
     */
    Scalar density(int phaseIdx) const
    {
        if (phaseIdx != phaseIndex)
            return 0;
        return density_;
    }

    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx != phaseIndex)
            return 0;
        return molarDensity_;
    }

    /*!
     * \brief Returns mean molar mass of a phase \f$\mathrm{[kg/mol]}\f$.
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     *
     * \param phaseIdx The index of the considered phase
     */
    Scalar meanMolarMass(int phaseIdx) const
    {
        if (phaseIdx != phaseIndex)
            return 0;
        return meanMolarMass_;
    };

    /*!
     * \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *
     * \param phaseIdx The index of the considered phase
     */
    Scalar phasePressure(int phaseIdx) const
    {
        return phasePressure_;
    }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

public:
    Scalar x1_;
    Scalar phasePressure_;
    Scalar density_;
    Scalar molarDensity_;
    Scalar meanMolarMass_;
    Scalar temperature_;
};

} // end namepace

#endif
