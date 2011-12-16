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
 * \brief Calcultes the fluid state from the primary variables in the
 *        Richards model.
 */
#ifndef DUMUX_RICHARDS_FLUID_STATE_HH
#define DUMUX_RICHARDS_FLUID_STATE_HH

#include <dumux/material/fluidstate.hh>

#include "richardsproperties.hh"

namespace Dumux
{
/*!
 * \brief Calcultes the fluid state from the primary variables in the
 *        Richards model.
 */
template <class TypeTag>
class RichardsFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                             RichardsFluidState<TypeTag> >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        pwIdx= Indices::pwIdx,
    };

public:
    //! The number of fluid phases used by the fluid state
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };

    /*!
     * \brief Updates the fluid quantities from the primary variables
     *        of the Richards model.
     *
     * \param pnRef The reference pressure of the non-wetting fluid phase \f$\mathrm{[Pa]}\f$
     * \param matParams The parameters for the capillary pressure/relative permeability material law
     * \param priVars The primary variables for which the fluid state ought to be calculated
     * \param temperature The temperature which should be used
     */
    void update(Scalar pnRef,
                const MaterialLawParams &matParams,
                const PrimaryVariables &priVars,
                Scalar temperature)
    {
        Scalar minPc = MaterialLaw::pC(matParams, 1.0);
        pressure_[wPhaseIdx] = priVars[pwIdx];
        pressure_[nPhaseIdx] = std::max(pnRef, priVars[pwIdx] + minPc);
        Sn_ = 1.0 - MaterialLaw::Sw(matParams, capillaryPressure());

        temperature_=temperature;
        
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(*this);
        densityWetting_ = FluidSystem::density(*this, paramCache, wPhaseIdx);
        viscosityWetting_ = FluidSystem::viscosity(*this, paramCache, wPhaseIdx);
    }

    /*!
     * \brief Returns the saturation [] of a phase.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
            return Scalar(1.0) - Sn_;
        else
            return Sn_;
    };

    /*!
     * \brief Returns the mass fraction [] of a component in a phase.
     *
     * \param phaseIdx The index of the fluid phase
     * \param compIdx The index of the chemical species
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        if (compIdx == phaseIdx)
            return 1.0;
        return 0;
    }

    /*!
     * \brief Returns the molar fraction [] of a component in a fluid phase.
     *
     * \param phaseIdx The index of the fluid phase
     * \param compIdx The index of the chemical species
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        return massFraction(phaseIdx, compIdx);
    }

    /*!
     * \brief Returns the molar density of a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param phaseIdx The index of the fluid phase
     * \param compIdx The index of the chemical species
     */
    Scalar molarDensity(int phaseIdx, int compIdx) const
    {
        return density(phaseIdx)/FluidSystem::molarMass(compIdx);
    };

    /*!
     * \brief Returns the molar density of a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param phaseIdx The index of the fluid phase
     * \param compIdx The index of the chemical species
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    {
        if (phaseIdx != compIdx)
            return 0;
        return molarDensity(phaseIdx);
    }

    /*!
     * \brief Returns the density of a phase \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar density(int phaseIdx) const
    { 
        if (phaseIdx == wPhaseIdx)
            return densityWetting_;
        return 1e-10;
    }

    /*!
     * \brief Returns the dynamic viscosity of a phase [Pa s].
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar viscosity(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
            return viscosityWetting_;
        return 1e-20; // something quite small
    }

    /*!
     * \brief Returns mean molar mass of a phase \f$\mathrm{[kg/mol]}\f$.
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return FluidSystem::molarMass(phaseIdx); };

    /*!
     * \brief Returns the partial pressure of a component in the gas phase \f$\mathrm{[Pa]}\f$.
     *
     * \param compIdx The index of the chemical species
     */
    Scalar fugacity(int compIdx) const
    {
        if (compIdx == wPhaseIdx)
            return 0;
        return pressure_[nPhaseIdx];
    }

    /*!
     * \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure \f$\mathrm{[Pa]}\f$
     */
    Scalar capillaryPressure() const
    { return pressure_[nPhaseIdx] - pressure_[wPhaseIdx]; }

    /*!
     * \brief Returns the temperature a single fluid \f$\mathrm{[K]}\f$.
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; };

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

public:
    Scalar pressure_[numPhases];
    Scalar densityWetting_;
    Scalar viscosityWetting_;
    Scalar temperature_;
    Scalar Sn_;
};

} // end namepace

#endif
