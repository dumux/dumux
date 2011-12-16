/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \brief Calculates the phase state from the primary variables in the
 *        2p2c model.
 */
#ifndef DUMUX_2P2C_FLUID_STATE_HH
#define DUMUX_2P2C_FLUID_STATE_HH

#include "2p2cproperties.hh"
#include "2p2cindices.hh"

#include <dumux/material/MpNcconstraintsolvers/computefromreferencephase.hh>
#include <dumux/material/MpNcconstraintsolvers/misciblemultiphasecomposition.hh>

#include <dumux/material/fluidstate.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief Calculates the phase state from the primary variables in the
 *        2p2c model.
 */
template <class TypeTag>
class TwoPTwoCFluidState
{
    typedef TwoPTwoCFluidState<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

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

    typedef Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

public:
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };

    /*!
     * \brief Update the phase state from the primary variables.
     *
     * \param primaryVars The primary variables
     * \param pcParams The parameters for the material law
     * \param temperature The temperature
     * \param phasePresence Stands either for nonwetting phase, wetting phase or both phases
     */
    template <class ParameterCache>
    void update(ParameterCache &paramCache,
                const PrimaryVariables &primaryVars,
                const MaterialLawParams &pcParams,
                int phasePresence)
    {
        Valgrind::CheckDefined(primaryVars);
        Valgrind::CheckDefined(temperature_);

        Valgrind::SetUndefined(moleFraction_);
        Valgrind::SetUndefined(density_);
        Valgrind::SetUndefined(avgMolarMass_);
        Valgrind::SetUndefined(pressure_);
        Valgrind::SetUndefined(Sg_);

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
            pressure_[lPhaseIdx] = primaryVars[pressureIdx];
            pressure_[gPhaseIdx] = pressure_[lPhaseIdx] + pC;
        }
        else if (formulation == pgSl) {
            pressure_[gPhaseIdx] = primaryVars[pressureIdx];
            pressure_[lPhaseIdx] = pressure_[gPhaseIdx] - pC;
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        Valgrind::CheckDefined(pressure_);

        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases) {
            // both phases are present, phase composition results from
            // the gas <-> liquid equilibrium. This is the job of the
            // "MiscibleMultiPhaseComposition" constraint solver
            MiscibleMultiPhaseComposition::solve(*this,
                                                 paramCache,
                                                 /*setViscosity=*/true,
                                                 /*setInternalEnergy=*/false);
            
        }
        else if (phasePresence == gPhaseOnly) {
            // only the gas phase is present, gas phase composition is
            // stored explicitly.

            // extract _mass_ (!) fractions in the gas phase, misuse
            // the moleFraction_ array for temporary storage
            moleFraction_[gPhaseIdx][lCompIdx] = primaryVars[switchIdx];
            moleFraction_[gPhaseIdx][gCompIdx] = 1 - moleFraction_[gPhaseIdx][lCompIdx];

            // calculate average molar mass of the gas phase
            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);
            Scalar X2 = moleFraction_[gPhaseIdx][gCompIdx]; // mass fraction of solvent in gas
            avgMolarMass_[gPhaseIdx] = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass to mole fractions
            moleFraction_[gPhaseIdx][lCompIdx] *= avgMolarMass_[gPhaseIdx]/M1;
            moleFraction_[gPhaseIdx][gCompIdx] *= avgMolarMass_[gPhaseIdx]/M2;

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(*this, 
                                             paramCache,
                                             gPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }
        else if (phasePresence == lPhaseOnly) {
            // only the gas phase is present, liquid phase composition is
            // stored explicitly.

            // extract _mass_ (!) fractions in the liquid phase, misuse
            // the moleFraction_ array for temporary storage
            moleFraction_[lPhaseIdx][gCompIdx] = primaryVars[switchIdx];
            moleFraction_[lPhaseIdx][lCompIdx] = 1 - moleFraction_[lPhaseIdx][gCompIdx];

            Scalar M1 = FluidSystem::molarMass(lCompIdx);
            Scalar M2 = FluidSystem::molarMass(gCompIdx);
            Scalar X2 = moleFraction_[lPhaseIdx][gCompIdx]; // mass fraction of solvent in gas
            avgMolarMass_[lPhaseIdx] = M1*M2/(M2 + X2*(M1 - M2));

            // convert mass to mole fractions
            moleFraction_[lPhaseIdx][lCompIdx] *= avgMolarMass_[lPhaseIdx]/M1;
            moleFraction_[lPhaseIdx][gCompIdx] *= avgMolarMass_[lPhaseIdx]/M2;

            // calculate the composition of the remaining phases. this
            // is the job of the "ComputeFromReferencePhase"
            // constraint solver
            ComputeFromReferencePhase::solve(*this, 
                                             paramCache,
                                             lPhaseIdx,
                                             /*setViscosity=*/true,
                                             /*setInternalEnergy=*/false);
        }

        Valgrind::CheckDefined(moleFraction_);
        Valgrind::CheckDefined(density_);
        Valgrind::CheckDefined(avgMolarMass_);
        Valgrind::CheckDefined(pressure_);
        Valgrind::CheckDefined(temperature_);
        Valgrind::CheckDefined(Sg_);
    }

    /*!
     * \brief Returns the saturation of a phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == lPhaseIdx)
            return Scalar(1.0) - Sg_;
        else
            return Sg_;
    };

    /*!
     * \brief Returns the molar density of a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    { return density_[phaseIdx]/avgMolarMass_[phaseIdx]; };

    /*!
     * \brief Returns the molar molarity of a component in a phase \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     *
     */
    Scalar molarity(int phaseIdx, int compIdx) const
    { return molarDensity(phaseIdx)*moleFraction(phaseIdx, compIdx); };

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     *
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        return
            moleFraction(phaseIdx, compIdx) 
            * FluidSystem::molarMass(compIdx)
            / avgMolarMass_[phaseIdx];
    }

    /*!
     * \brief Returns the mole fraction of a component in a fluid phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        return moleFraction_[phaseIdx][compIdx];
    }

    /*!
     * \brief Returns the mass density of a phase \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief Returns the dynamic viscosity of a phase \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*!
     * \brief Returns the specific internal energy of a phase \f$\mathrm{[J/kg]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(int phaseIdx) const
    { return internalEnergy_[phaseIdx]; }

    /*!
     * \brief Returns the specific enthalpy of a phase \f$\mathrm{[J/kg]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(int phaseIdx) const
    { return internalEnergy_[phaseIdx] + pressure_[phaseIdx]/density_[phaseIdx]; }

    /*!
     * \brief Returns mean molar mass of a phase \f$\mathrm{[kg/mol]}\f$.
     *
     * \param phaseIdx The phase index
     *
     * This is equivalent to the sum of all component molar masses
     * weighted by their respective mole fraction.
     */
    Scalar averageMolarMass(int phaseIdx) const
    { return avgMolarMass_[phaseIdx]; };

    /*!
     * \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return pressure_[phaseIdx]; }

    /*!
     * \brief The fugacity of a component in a phase [Pa]
     */
    Scalar fugacity(int phaseIdx, int compIdx) const
    { return fugacityCoefficient(phaseIdx, compIdx)*moleFraction(phaseIdx, compIdx)*pressure(phaseIdx); }

    /*!
     * \brief The fugacity coefficient of a component in a phase [Pa]
     */
    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    { return fugacityCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief Returns the capillary pressure \f$\mathrm{[Pa]}\f$
     */
    Scalar capillaryPressure() const
    { return pressure_[gPhaseIdx] - pressure_[lPhaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

    /*!
     * \brief Returns the temperature of a fluid phases \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; };

    /*!
     * \brief Set the mole fraction of a component in a phase []
     */
    void setMoleFraction(int phaseIdx, int compIdx, Scalar value)
    { moleFraction_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Set the fugacity of a component in a phase []
     */
    void setFugacityCoefficient(int phaseIdx, int compIdx, Scalar value)
    { fugacityCoefficient_[phaseIdx][compIdx] = value; }

    /*!
     * \brief Set the density of a phase [kg / m^3]
     */
    void setDensity(int phaseIdx, Scalar value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Set the temperature of *all* phases [K]
     */
    void setTemperature(Scalar value)
    { temperature_ = value; }

    /*!
     * \brief Set the specific internal energy of a phase [J/m^3]
     */
    void setInternalEnergy(int phaseIdx, Scalar value)
    { internalEnergy_[phaseIdx] = value; }

    /*!
     * \brief Set the dynamic viscosity of a phase [Pa s]
     */
    void setViscosity(int phaseIdx, Scalar value)
    { viscosity_[phaseIdx] = value; }

    /*!
     * \brief Calculatate the mean molar mass of a phase given that
     *        all mole fractions have been set
     */
    void updateAverageMolarMass(int phaseIdx)
    {
        avgMolarMass_[phaseIdx] = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            avgMolarMass_[phaseIdx] += moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);
        }
        Valgrind::CheckDefined(avgMolarMass_[phaseIdx]);
    }

public:
    Scalar moleFraction_[numPhases][numComponents];
    Scalar fugacityCoefficient_[numPhases][numComponents];
    Scalar internalEnergy_[numPhases];
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];
    Scalar avgMolarMass_[numPhases];
    Scalar pressure_[numPhases];
    Scalar temperature_;
    Scalar Sg_;
};

} // end namepace

#endif
