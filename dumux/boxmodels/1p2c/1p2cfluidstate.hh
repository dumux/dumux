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

namespace Dumux
{
/*!
 * \brief Calcultes the phase state from the primary variables in the
 *        1p2c model.
 */
template <class TypeTag>
class OnePTwoCFluidState
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

        phaseIdx = Indices::phaseIdx,
        comp0Idx = Indices::comp0Idx,
        comp1Idx = Indices::comp1Idx,
    };

    static const bool useMoles = GET_PROP_VALUE(TypeTag, PTAG(UseMoles));

public:
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };

    /*!
     * \brief Update the phase state from the primary variables.
     *
     * \param primaryVars The primary variables
     */
    template <class ParameterCache>
    void update(ParameterCache &paramCache,
                const PrimaryVariables &primaryVars)
    {
        Valgrind::CheckDefined(primaryVars);
        Valgrind::CheckDefined(temperature_);

        pressure_ =  primaryVars[pressureIdx];
        x1_ = primaryVars[x1Idx]; //mole or mass fraction of component 1

        if(!useMoles) //mass-fraction formulation
        {
            Scalar M0 = FluidSystem::molarMass(comp0Idx);
            Scalar M1 = FluidSystem::molarMass(comp1Idx);
            //meanMolarMass if x1_ is a massfraction
            meanMolarMass_ = M0*M1/(M1 + x1_*(M0 - M1));
        }
        else //mole-fraction formulation
        {
            //meanMolarMass if x1_ is a molefraction
            meanMolarMass_ =
            (1 - x1_)*FluidSystem::molarMass(comp0Idx) +
            (x1_    )*FluidSystem::molarMass(comp1Idx);
        }

        paramCache.updatePhase(*this, /*phaseIdx=*/0);
        
        density_ = FluidSystem::density(*this, paramCache, /*phaseIdx=*/0);
        viscosity_ = FluidSystem::viscosity(*this, paramCache, /*phaseIdx=*/0);

        Valgrind::CheckDefined(x1_);
        Valgrind::CheckDefined(pressure_);
        Valgrind::CheckDefined(density_);
        Valgrind::CheckDefined(viscosity_);
        Valgrind::CheckDefined(meanMolarMass_);
        Valgrind::CheckDefined(temperature_);
    }

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     *
     * \param phaseIndex The index of the considered phase
     * \param compIdx The index of the considered component
     */
    Scalar moleFraction(int phaseIndex, int compIdx) const
    {
        if(!useMoles) //mass-fraction formulation
        {
            ///if x1_ is a massfraction
            Scalar moleFrac1(x1_);
            moleFrac1 *= meanMolarMass_/FluidSystem::molarMass(comp1Idx);

            if (compIdx==comp0Idx)
                return 1-moleFrac1;
            else if (compIdx==comp1Idx)
                return moleFrac1;
            else
                return 0.0;
        }
        else //mole-fraction formulation
        {
            //if x1_ is a molefraction
            if (compIdx==comp0Idx)
                return 1-x1_;
            else if (compIdx==comp1Idx)
                return x1_;
            else
                return 0.0;
        }
    }

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     *
     * \param phaseIndex The index of the considered phase
     * \param compIdx The index of the considered component
     */
    Scalar massFraction(int phaseIndex, int compIdx) const
    {
        if(!useMoles)
        {
            //if x1_ is a mass fraction
            if (compIdx==comp0Idx)
                return 1-x1_;
            else if (compIdx==comp1Idx)
                return x1_;
            return 0.0;
        }
        else
        {
            //if x1_ is a molefraction
            return
                moleFraction(phaseIndex, compIdx)*
                FluidSystem::molarMass(compIdx)
                / meanMolarMass_;
        }
    }

    /*!
     * \brief Returns the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIndex The index of the considered phase
     *
     */
    Scalar density(int phaseIndex) const
    {
        return density_;
    }

    /*!
     * \brief Returns the dynamic viscosity of a phase \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_; }

    /*!
     * \brief Returns the molar density of a phase \f$\mathrm{[mole/m^3]}\f$.
     *
     * \param phaseIndex The index of the considered phase
     *
     */
    Scalar molarDensity(int phaseIndex) const
    {
        return density_/meanMolarMass_;
    }

    /*!
     * \brief Returns the molar concentration of a component in a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param phaseIndex The index of the considered phase
     * \param compIdx The index of the considered component
     */
    Scalar molarity(int phaseIndex, int compIdx) const
    { return molarDensity(phaseIndex)*moleFraction(phaseIndex, compIdx); };

    /*!
     * \brief Returns mean molar mass of a phase \f$\mathrm{[kg/mol]}\f$.
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     *
     * \param phaseIndex The index of the considered phase
     */
    Scalar averageMolarMass(int phaseIndex) const
    {
        return meanMolarMass_;
    };

    /*!
     * \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *
     * \param phaseIndex The index of the considered phase
     */
    Scalar pressure(int phaseIndex) const
    {
        return pressure_;
    }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

    /*!
     * \brief Set the temperature of *all* phases [K]
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

public:
    Scalar x1_;
    Scalar pressure_;
    Scalar density_;
    Scalar viscosity_;
    Scalar meanMolarMass_;
    Scalar temperature_;
};

} // end namepace

#endif
