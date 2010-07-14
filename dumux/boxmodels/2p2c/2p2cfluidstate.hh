// $Id: 2p2cfluidstate.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Calcultes the phase state from the primary variables in the
 *        2pNc model.
 */
#ifndef DUMUX_2P2C_PHASE_STATE_HH
#define DUMUX_2P2C_PHASE_STATE_HH

#include "2p2cproperties.hh"

#include <dumux/material/fluidstate.hh>

namespace Dumux
{
/*!
 * \brief Calcultes the phase state from the primary variables in the
 *        2pNc model.
 */
template <class TypeTag>
class TwoPTwoCFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                             TwoPTwoCFluidState<TypeTag> >
{
    typedef TwoPTwoCFluidState<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename MaterialLaw::Params                       MaterialLawParams;

    enum {
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,

        switchIdx = Indices::switchIdx,
        pressureIdx = Indices::pressureIdx,
    };

    // present phases
    enum {
        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    // formulations
    enum {
        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)),
        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl
    };

public:
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };
    enum { numSolvents = 1 };

    /*!
     * \brief Update the phase state from the primary variables.
     */
    void update(const PrimaryVarVector &primaryVars,
                const MaterialLawParams &pcParams,
                Scalar temperature,
                int phasePresence)
    {
        Valgrind::CheckDefined(primaryVars);

        Valgrind::SetUndefined(concentration_);
        Valgrind::SetUndefined(totalConcentration_);
        Valgrind::SetUndefined(avgMolarMass_);
        Valgrind::SetUndefined(phasePressure_);
        Valgrind::SetUndefined(temperature_);
        Valgrind::SetUndefined(Sg_);

        temperature_ = temperature;
        Valgrind::CheckDefined(temperature_);

        // extract non-wetting phase pressure
        if (phasePresence == gPhaseOnly)
            Sg_ = 1.0;
        else if (phasePresence == lPhaseOnly) {
            Sg_ = 0.0;
        }
        else if (phasePresence == bothPhases) {
            if (formulation == plSg)
                Sg_ = primaryVars[switchIdx];
            else if (formulation == pgSl)
                Sg_ = 1.0 - primaryVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        Valgrind::CheckDefined(Sg_);

        // calculate capillary pressure
        Scalar pC = MaterialLaw::pC(pcParams, 1 - Sg_);

        // extract the pressures
        if (formulation == plSg) {
            phasePressure_[lPhaseIdx] = primaryVars[pressureIdx];
            phasePressure_[gPhaseIdx] = phasePressure_[lPhaseIdx] + pC;
        }
        else if (formulation == pgSl) {
            phasePressure_[gPhaseIdx] = primaryVars[pressureIdx];
            phasePressure_[lPhaseIdx] = phasePressure_[gPhaseIdx] - pC;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        Valgrind::CheckDefined(phasePressure_);

        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases) {
            // both phases are present, phase composition results from
            // the gas <-> liquid equilibrium.
            FluidSystem::computeEquilibrium(*this, -1);
        }
        else if (phasePresence == gPhaseOnly) {
            // only the gas phase is present, gas phase composition is
            // stored explicitly.

            // extract _mass_ (!) fractions in the gas phase, misuse
            // the concentration_ array for temporary storage
            concentration_[gPhaseIdx][lCompIdx] = primaryVars[switchIdx];
            concentration_[gPhaseIdx][gCompIdx] = 1 - concentration_[gPhaseIdx][lCompIdx];

            // calculate average molar mass of the gas phase
            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);
            Scalar X2 = concentration_[gPhaseIdx][gCompIdx]; // mass fraction of solvent in gas
            avgMolarMass_[gPhaseIdx] = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass to mole fractions
            concentration_[gPhaseIdx][lCompIdx] *= avgMolarMass_[gPhaseIdx]/M1;
            concentration_[gPhaseIdx][gCompIdx] *= avgMolarMass_[gPhaseIdx]/M2;

            // convert to real concentrations
            totalConcentration_[gPhaseIdx] = 1.0;
            Scalar rhog = FluidSystem::phaseDensity(gPhaseIdx,
                                                    temperature_,
                                                    phasePressure_[gPhaseIdx],
                                                    *this);
            totalConcentration_[gPhaseIdx] = rhog / avgMolarMass_[gPhaseIdx];
            concentration_[gPhaseIdx][lCompIdx] *= totalConcentration_[gPhaseIdx];
            concentration_[gPhaseIdx][gCompIdx] *= totalConcentration_[gPhaseIdx];

            // tell the fluid system to calculate the composition of
            // the remaining phases
            FluidSystem::computeEquilibrium(*this, gPhaseIdx);
        }
        else if (phasePresence == lPhaseOnly) {
            // only the gas phase is present, liquid phase composition is
            // stored explicitly.

            // extract _mass_ (!) fractions in the liquid phase, misuse
            // the concentration_ array for temporary storage
            concentration_[lPhaseIdx][gCompIdx] = primaryVars[switchIdx];
            concentration_[lPhaseIdx][lCompIdx] = 1 - concentration_[lPhaseIdx][gCompIdx];

            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);
            Scalar X2 = concentration_[lPhaseIdx][gCompIdx]; // mass fraction of solvent in gas
            avgMolarMass_[lPhaseIdx] = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass to mole fractions
            concentration_[lPhaseIdx][lCompIdx] *= avgMolarMass_[lPhaseIdx]/M1;
            concentration_[lPhaseIdx][gCompIdx] *= avgMolarMass_[lPhaseIdx]/M2;

            // convert to real concentrations
            totalConcentration_[lPhaseIdx] = 1.0;
            Scalar rhol = FluidSystem::phaseDensity(lPhaseIdx,
                                                    temperature_,
                                                    phasePressure_[lPhaseIdx],
                                                    *this);
            totalConcentration_[lPhaseIdx] = rhol / avgMolarMass_[lPhaseIdx];
            concentration_[lPhaseIdx][lCompIdx] *= totalConcentration_[lPhaseIdx];
            concentration_[lPhaseIdx][gCompIdx] *= totalConcentration_[lPhaseIdx];

            // tell the fluid system to calculate the composition of
            // the remaining phases
            FluidSystem::computeEquilibrium(*this, lPhaseIdx);
        }

        Valgrind::CheckDefined(concentration_);
        Valgrind::CheckDefined(totalConcentration_);
        Valgrind::CheckDefined(avgMolarMass_);
        Valgrind::CheckDefined(phasePressure_);
        Valgrind::CheckDefined(temperature_);
        Valgrind::CheckDefined(Sg_);
    }

public:
    /*!
     * \brief Retrieves the phase composition and pressure from a
     *        phase composition class.
     *
     * This method is called by the fluid system's
     * computeEquilibrium()
     */
    template <class PhaseCompo>
    void assignPhase(int phaseIdx, const PhaseCompo &compo)
    {
        avgMolarMass_[phaseIdx] = compo.meanMolarMass();
        phasePressure_[phaseIdx] = compo.pressure();

        totalConcentration_[phaseIdx] = compo.totalConcentration();
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            concentration_[phaseIdx][compIdx] = compo.concentration(compIdx);
    };

    /*!
     * \brief Returns the saturation of a phase.
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == lPhaseIdx)
            return Scalar(1.0) - Sg_;
        else
            return Sg_;
    };

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    {
        return
            concentration_[phaseIdx][compIdx]/
            totalConcentration_[phaseIdx];
    }

    /*!
     * \brief Returns the total concentration of a phase [mol / m^3].
     *
     * This is equivalent to the sum of all component concentrations.
     */
    Scalar totalConcentration(int phaseIdx) const
    { return totalConcentration_[phaseIdx]; };

    /*!
     * \brief Returns the concentration of a component in a phase [mol / m^3].
     */
    Scalar concentration(int phaseIdx, int compIdx) const
    { return concentration_[phaseIdx][compIdx]; };

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    {
        return
            moleFrac(phaseIdx, compIdx) *
            FluidSystem::molarMass(compIdx)/avgMolarMass_[phaseIdx];
    }

    /*!
     * \brief Returns the density of a phase [kg / m^3].
     */
    Scalar density(int phaseIdx) const
    { return totalConcentration_[phaseIdx]*avgMolarMass_[phaseIdx]; }

    /*!
     * \brief Returns mean molar mass of a phase [kg / mol].
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return avgMolarMass_[phaseIdx]; };

    /*!
     * \brief Returns the pressure of a fluid phase [Pa].
     */
    Scalar phasePressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the capillary pressure [Pa]
     */
    Scalar capillaryPressure() const
    { return phasePressure_[gPhaseIdx] - phasePressure_[lPhaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids [K].
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

public:
    Scalar concentration_[numPhases][numComponents];
    Scalar totalConcentration_[numPhases];
    Scalar avgMolarMass_[numPhases];
    Scalar phasePressure_[numPhases];
    Scalar temperature_;
    Scalar Sg_;
};

} // end namepace

#endif
