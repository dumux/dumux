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
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */
#ifndef DUMUX_2P2C_VOLUME_VARIABLES_HH
#define DUMUX_2P2C_VOLUME_VARIABLES_HH

#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */
template <class Traits>
class TwoPTwoCVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables< Traits >;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPTwoCVolumeVariables<Traits> >;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;

    static constexpr int numFluidComps = ParentType::numComponents();
    // component indices
    enum
    {
        comp0Idx = Traits::FluidSystem::comp0Idx,
        comp1Idx = Traits::FluidSystem::comp1Idx,
        phase0Idx = Traits::FluidSystem::phase0Idx,
        phase1Idx = Traits::FluidSystem::phase1Idx
    };

    // phase presence indices
    enum
    {
        firstPhaseOnly = ModelTraits::Indices::firstPhaseOnly,
        secondPhaseOnly = ModelTraits::Indices::secondPhaseOnly,
        bothPhases = ModelTraits::Indices::bothPhases
    };

    // primary variable indices
    enum
    {
        switchIdx = ModelTraits::Indices::switchIdx,
        pressureIdx = ModelTraits::Indices::pressureIdx
    };

    // formulations
    static constexpr auto formulation = ModelTraits::priVarFormulation();

    // further specifications on the variables update
    static constexpr bool useConstraintSolver = Traits::useConstraintSolver;

    using PermeabilityType = typename Traits::PermeabilityType;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, typename Traits::FluidSystem>;
    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition< Scalar, typename Traits::FluidSystem >;
public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename Traits::FluidState;
    //! The fluid system used here
    using FluidSystem = typename Traits::FluidSystem;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    //! return whether moles or masses are balanced
    static constexpr bool useMoles() { return ModelTraits::useMoles(); }
    //! return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }

    // check for permissive combinations
    static_assert(useMoles() || (!useMoles() && useConstraintSolver), "if !UseMoles, UseConstraintSolver has to be set to true");
    static_assert(ModelTraits::numPhases() == 2, "NumPhases set in the model is not two!");
    static_assert(ModelTraits::numComponents() == 2, "NumComponents set in the model is not two!");
    static_assert((formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0), "Chosen TwoPFormulation not supported!");

    // The computations in the explicit composition update most probably assume a liquid-gas interface with
    // liquid as first phase. TODO: is this really needed? The constraint solver does the job anyway, doesn't it?
    static_assert(useConstraintSolver || (!FluidSystem::isGas(phase0Idx) && FluidSystem::isGas(phase1Idx)),
                   "Explicit composition calculation has to be re-checked for NON-liquid-gas equilibria");

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol, const Problem& problem, const Element& element, const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        // Second instance of a parameter cache. Could be avoided if
        // diffusion coefficients also became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& matParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        const int nPhaseIdx = 1 - wPhaseIdx;

        // relative permeabilities -> require wetting phase saturation as parameter!
        relativePermeability_[wPhaseIdx] = MaterialLaw::krw(matParams, saturation(wPhaseIdx));
        relativePermeability_[nPhaseIdx] = MaterialLaw::krn(matParams, saturation(wPhaseIdx));

        // binary diffusion coefficients
        diffCoeff_[phase0Idx] = FluidSystem::binaryDiffusionCoefficient(fluidState_, paramCache, phase0Idx, comp0Idx, comp1Idx);
        diffCoeff_[phase1Idx] = FluidSystem::binaryDiffusionCoefficient(fluidState_, paramCache, phase1Idx, comp0Idx, comp1Idx);

        // porosity & permeabilty
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
    }

    /*!
     * \brief Complete the fluid state
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The problem
     * \param element The element
     * \param scv The sub control volume
     * \param fluidState The fluid state
     *
     * Set temperature, saturations, capillary pressures, viscosities, densities and enthalpies.
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);

        // set the saturations
        if (phasePresence == firstPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 1.0);
            fluidState.setSaturation(phase1Idx, 0.0);
        }
        else if (phasePresence == secondPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 0.0);
            fluidState.setSaturation(phase1Idx, 1.0);
        }
        else if (phasePresence == bothPhases)
        {
            if (formulation == TwoPFormulation::p0s1)
            {
                fluidState.setSaturation(phase1Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase0Idx, 1 - priVars[switchIdx]);
            }
            else
            {
                fluidState.setSaturation(phase0Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase1Idx, 1 - priVars[switchIdx]);
            }
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase presence.");

        // set pressures of the fluid phases
        pc_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx));
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase1Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] + pc_
                                                                       : priVars[pressureIdx] - pc_);
        }
        else
        {
            fluidState.setPressure(phase1Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase0Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] - pc_
                                                                       : priVars[pressureIdx] + pc_);
        }

        // calculate the phase compositions
        typename FluidSystem::ParameterCache paramCache;

        // If constraint solver is not used, get the phase pressures and set the fugacity coefficients here
        if(!useConstraintSolver)
        {
            for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++ phaseIdx)
            {
                assert(FluidSystem::isIdealMixture(phaseIdx));
                for (int compIdx = 0; compIdx < ModelTraits::numComponents(); ++ compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }

        // now comes the tricky part: calculate phase compositions
        const Scalar p0 = fluidState.pressure(phase0Idx);
        const Scalar p1 = fluidState.pressure(phase1Idx);
        if (phasePresence == bothPhases)
        {
            // both phases are present, phase compositions are a result
            // of the equilibrium between the phases. This is the job
            // of the "MiscibleMultiPhaseComposition" constraint solver
            if(useConstraintSolver)
                MiscibleMultiPhaseComposition::solve(fluidState,
                                                     paramCache);
            // ... or calculated explicitly this way ...
            else
            {
                // get the partial pressure of the main component of the first phase within the
                // second phase == vapor pressure due to equilibrium. Note that in this case the
                // fugacityCoefficient * p is the vapor pressure (see implementation in respective fluidsystem)
                const Scalar partPressLiquid = FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx)*p0;

                // get the partial pressure of the main component of the gas phase
                const Scalar partPressGas = p1 - partPressLiquid;

                // calculate the mole fractions of the components within the nonwetting phase
                const Scalar xnn = partPressGas / p1;
                const Scalar xnw = partPressLiquid / p1;

                // calculate the mole fractions of the components within the wetting phase
                // note that in this case the fugacityCoefficient * p is the Henry Coefficient
                // (see implementation in respective fluidsystem)
                const Scalar xwn = partPressGas / (FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx)*p0);
                const Scalar xww = 1.0 - xwn;

                // set all mole fractions
                fluidState.setMoleFraction(phase0Idx, comp0Idx, xww);
                fluidState.setMoleFraction(phase0Idx, comp1Idx, xwn);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, xnw);
                fluidState.setMoleFraction(phase1Idx, comp1Idx, xnn);
            }
        }
        else if (phasePresence == secondPhaseOnly)
        {
            // only the second phase is present, composition is stored explicitly.
            if( useMoles() )
            {
                fluidState.setMoleFraction(phase1Idx, comp1Idx, 1 - priVars[switchIdx]);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, priVars[switchIdx]);
            }
            // setMassFraction() has only to be called 1-numComponents times
            else
                fluidState.setMassFraction(phase1Idx, comp0Idx, priVars[switchIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver)
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 phase1Idx);
            // ... or calculated explicitly this way ...
            else
            {
                // note that the water phase is actually not existing!
                // thus, this is used as phase switch criterion
                const Scalar xnw = priVars[switchIdx];
                const Scalar xnn = 1.0 - xnw;

                // first, xww:
                // xnw * pn = "actual" (hypothetical) vapor pressure
                // fugacityCoefficient * pw = vapor pressure given by thermodynamic conditions
                // Here, xww is not actually the mole fraction of water in the wetting phase
                // xww is only the ratio of "actual" vapor pressure / "thermodynamic" vapor pressure
                // If xww > 1 : gas is over-saturated with water vapor,
                // condensation takes place (see switch criterion in model)
                const Scalar xww = xnw*p1/( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx)*p0 );

                // second, xwn:
                // partialPressure / xwn = Henry
                // partialPressure = xnn * pn
                // xwn = xnn * pn / Henry
                // Henry = fugacityCoefficient * pw
                const Scalar xwn = xnn*p1/( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx)*p0 );

                fluidState.setMoleFraction(phase0Idx, comp0Idx, xww);
                fluidState.setMoleFraction(phase0Idx, comp1Idx, xwn);
            }
        }
        else if (phasePresence == firstPhaseOnly)
        {
            // only the wetting phase is present, i.e. wetting phase
            // composition is stored explicitly.
            if( useMoles() ) // mole-fraction formulation
            {
                fluidState.setMoleFraction(phase0Idx, comp0Idx, 1-priVars[switchIdx]);
                fluidState.setMoleFraction(phase0Idx, comp1Idx, priVars[switchIdx]);
            }
            // setMassFraction() has only to be called 1-numComponents times
            else // mass-fraction formulation
                fluidState.setMassFraction(phase0Idx, comp1Idx, priVars[switchIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver)
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 phase0Idx);
            // ... or calculated explicitly this way ...
            else
            {
                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                const Scalar xwn = priVars[switchIdx];

                // first, xnw:
                // psteam = xnw * pn = partial pressure of water in gas phase
                // psteam = fugacityCoefficient * pw
                const Scalar xnw = ( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx)*p0 )/p1;

                // second, xnn:
                // xwn = partialPressure / Henry
                // partialPressure = pn * xnn
                // xwn = pn * xnn / Henry
                // xnn = xwn * Henry / pn
                // Henry = fugacityCoefficient * pw
                const Scalar xnn = xwn*( FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx)*p0 )/p1;

                fluidState.setMoleFraction(phase1Idx, comp1Idx, xnn);
                fluidState.setMoleFraction(phase1Idx, comp0Idx, xnw);
            }
        }

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            // set the viscosity and desity here if constraintsolver is not used
            if(!useConstraintSolver)
            {
                paramCache.updateComposition(fluidState, phaseIdx);
                const Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);
                Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
                fluidState.setMolarDensity(phaseIdx, rhoMolar);
            }

            // compute and set the enthalpy
            const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx,mu);
            Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }

    /*!
     * \brief Returns the phase state within the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(const int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar massFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.massFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mole fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    { return fluidState_.moleFraction(phaseIdx, compIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[kg/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity of the fluid within the
     *        control volume in \f$\mathrm{[Pa s]}\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(const int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[mol/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx); }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature within the control volume in \f$[K]\f$.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the relative permeability of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar relativePermeability(const int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    { return relativePermeability_[phaseIdx]/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return pc_; }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the average permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        if(phaseIdx == compIdx)
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient called for phaseIdx = compIdx");
        else
            return diffCoeff_[phaseIdx];
    }

private:
    FluidState fluidState_;
    SolidState solidState_;

    Scalar pc_;                     //!< The capillary pressure
    PermeabilityType permeability_; //!< Effective permeability within the control volume

    //!< Relative permeability within the control volume
    std::array<Scalar, ModelTraits::numPhases()> relativePermeability_;

    //!< Binary diffusion coefficients for the phases
    std::array<Scalar, ModelTraits::numPhases()> diffCoeff_;
};

} // end namespace Dumux

#endif
