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
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/discretization/methods.hh>

#include "indices.hh" // for formulation

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase two-component model.
 */
template <class TypeTag>
class TwoPTwoCVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    // component indices
    enum
    {
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    // present phases
    enum
    {
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    // formulations
    enum
    {
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        pwsn = TwoPTwoCFormulation::pwsn,
        pnsw = TwoPTwoCFormulation::pnsw
    };

    // primary variable indices
    enum
    {
        switchIdx = Indices::switchIdx,
        pressureIdx = Indices::pressureIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FluidSystem>;

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr bool useKelvinEquation = GET_PROP_VALUE(TypeTag, UseKelvinEquation);
    static constexpr bool useConstraintSolver = GET_PROP_VALUE(TypeTag, UseConstraintSolver);
    static_assert(useMoles || (!useMoles && useConstraintSolver),
                  "if UseMoles is set false, UseConstraintSolver has to be set to true");

    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem, useKelvinEquation>;

    static constexpr int numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases();
    static constexpr int numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents();
public:
    //! The type of the object returned by the fluidState() method
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    //! Pull member functions of the parent for deriving classes
    using ParentType::temperature;
    using ParentType::enthalpy;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        Implementation::completeFluidState(elemSol, problem, element, scv, fluidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // relative permeabilities
            Scalar kr;
            if (phaseIdx == wPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(wPhaseIdx));
            else // ATTENTION: krn requires the wetting phase saturation
                // as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));
            relativePermeability_[phaseIdx] = kr;
            Valgrind::CheckDefined(relativePermeability_[phaseIdx]);

            // binary diffusion coefficients
            diffCoeff_[phaseIdx] = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                           paramCache,
                                                                           phaseIdx,
                                                                           wCompIdx,
                                                                           nCompIdx);
        }

        // porosity & permeabilty
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
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
    template<class ElementSolution>
    static void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);
        const auto phasePresence = priVars.state();

        /////////////
        // set the saturations
        /////////////
        Scalar sn;
        if (phasePresence == nPhaseOnly)
            sn = 1.0;
        else if (phasePresence == wPhaseOnly) {
            sn = 0.0;
        }
        else if (phasePresence == bothPhases) {
            if (formulation == pwsn)
                sn = priVars[switchIdx];
            else if (formulation == pnsw)
                sn = 1.0 - priVars[switchIdx];
            else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        fluidState.setSaturation(wPhaseIdx, 1 - sn);
        fluidState.setSaturation(nPhaseIdx, sn);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        Scalar pc = MaterialLaw::pc(materialParams, 1 - sn);

        if (formulation == pwsn) {
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx] + pc);
        }
        else if (formulation == pnsw) {
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx]);
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx] - pc);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");

        /////////////
        // calculate the phase compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;

        //get the phase pressures and set the fugacity coefficients here if constraintsolver is not used
        Scalar pn = 0;
        Scalar pw = 0;

        if(!useConstraintSolver) {
            if (formulation == pwsn) {
                pw = priVars[pressureIdx];
                pn = pw + pc;
            }
            else {
                pn = priVars[pressureIdx];
                pw = pn - pc;
            }

            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                assert(FluidSystem::isIdealMixture(phaseIdx));

                for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }

        // now comes the tricky part: calculate phase compositions
        if (phasePresence == bothPhases) {
            // both phases are present, phase compositions are a
            // result of the nonwetting <-> wetting equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver
            if(useConstraintSolver) {
                MiscibleMultiPhaseComposition::solve(fluidState,
                                                     paramCache,
                                                     /*setViscosity=*/true,
                                                     /*setEnthalpy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                //get the partial pressure of the main component of the the wetting phase ("H20") within the nonwetting (gas) phase == vapor pressure due to equilibrium
                //note that in this case the fugacityCoefficient * pw is the vapor pressure (see implementation in respective fluidsystem)
                Scalar partPressH2O = FluidSystem::fugacityCoefficient(fluidState,
                                                                       wPhaseIdx,
                                                                       wCompIdx) * pw;

                if (useKelvinEquation)
                    partPressH2O = FluidSystem::kelvinVaporPressure(fluidState, wPhaseIdx, wCompIdx);

                // get the partial pressure of the main component of the the nonwetting (gas) phase ("Air")
                Scalar partPressAir = pn - partPressH2O;

                //calculate the mole fractions of the components within the nonwetting phase
                Scalar xnn = partPressAir/pn;
                Scalar xnw = partPressH2O/pn;

                // calculate the mole fractions of the components within the wetting phase
                //note that in this case the fugacityCoefficient * pw is the Henry Coefficient (see implementation in respective fluidsystem)
                Scalar xwn = partPressAir
                  / (FluidSystem::fugacityCoefficient(fluidState,
                                                      wPhaseIdx,nCompIdx)
                  * pw);

                Scalar xww = 1.0 -xwn;

                //set all mole fractions
                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, xnn);
            }
        }
        else if (phasePresence == nPhaseOnly)
        {
            // only the nonwetting phase is present, i.e. nonwetting phase
            // composition is stored explicitly.
            if(useMoles)
            {
                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1 - priVars[switchIdx]);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, priVars[switchIdx]);
            }
            else
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(nPhaseIdx, wCompIdx, priVars[switchIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver) {
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 nPhaseIdx,
                                                 /*setViscosity=*/true,
                                                 /*setEnthalpy=*/false);
            }
            // ... or calculated explicitly this way ...
            else {
                // note that the water phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xnw = priVars[switchIdx];
                Scalar xnn = 1.0 -xnw;

                //first, xww:
                // xnw * pn = "actual" (hypothetical) vapor pressure
                // fugacityCoefficient * pw = vapor pressure given by thermodynamic conditions
                // Here, xww is not actually the mole fraction of water in the wetting phase
                // xww is only the ratio of "actual" vapor pressure / "thermodynamic" vapor pressure
                // If xww > 1 : gas is over-saturated with water vapor,
                // condensation takes place (see switch criterion in model)
                Scalar xww = xnw * pn
                  / (FluidSystem::fugacityCoefficient(fluidState,
                                                      wPhaseIdx,wCompIdx)
                     * pw);

                // now, xwn:
                //partialPressure / xwn = Henry
                //partialPressure = xnn * pn
                //xwn = xnn * pn / Henry
                // Henry = fugacityCoefficient * pw
                Scalar xwn = xnn * pn / (FluidSystem::fugacityCoefficient(fluidState,
                                                                          wPhaseIdx,nCompIdx)
                                         * pw);

                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xww);
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            // only the wetting phase is present, i.e. wetting phase
            // composition is stored explicitly.
            if(useMoles) // mole-fraction formulation
            {
                fluidState.setMoleFraction(wPhaseIdx, wCompIdx, 1-priVars[switchIdx]);
                fluidState.setMoleFraction(wPhaseIdx, nCompIdx, priVars[switchIdx]);
            }
            else // mass-fraction formulation
            {
                // setMassFraction() has only to be called 1-numComponents times
                fluidState.setMassFraction(wPhaseIdx, nCompIdx, priVars[switchIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). This is the job
            // of the "ComputeFromReferencePhase" constraint solver
            if (useConstraintSolver) {
                ComputeFromReferencePhase::solve(fluidState,
                                                 paramCache,
                                                 wPhaseIdx,
                                                 /*setViscosity=*/true,
                                                 /*setEnthalpy=*/false);
            }
            // ... or calculated explicitly this way ...
            else
            {
                // note that the gas phase is actually not existing!
                // thus, this is used as phase switch criterion
                Scalar xwn = priVars[switchIdx];

                //first, xnw:
                //psteam = xnw * pn = partial pressure of water in gas phase
                //psteam = fugacityCoefficient * pw
                Scalar xnw = (FluidSystem::fugacityCoefficient(fluidState,
                                                               wPhaseIdx,wCompIdx)
                              * pw) / pn ;

                //now, xnn:
                // xwn = partialPressure / Henry
                // partialPressure = pn * xnn
                // xwn = pn * xnn / Henry
                // xnn = xwn * Henry / pn
                // Henry = fugacityCoefficient * pw
                Scalar xnn = xwn * (FluidSystem::fugacityCoefficient(fluidState,
                                                                     wPhaseIdx,nCompIdx)
                                    * pw) / pn ;

                fluidState.setMoleFraction(nPhaseIdx, nCompIdx, xnn);
                fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
            }
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // set the viscosity and desity here if constraintsolver is not used
            if(!useConstraintSolver)
            {
                paramCache.updateComposition(fluidState, phaseIdx);
                const Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);
                const Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
                fluidState.setViscosity(phaseIdx,mu);
            }
            // compute and set the enthalpy
            Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }


    /*!
     * \brief Returns the phase state within the control volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

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
    { return fluidState_.density(phaseIdx) / fluidState_.averageMolarMass(phaseIdx); }

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
    {
        return relativePermeability_[phaseIdx];
    }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume in \f$[s*m/kg]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(const int phaseIdx) const
    {
        return relativePermeability_[phaseIdx]/fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume in \f$[-]\f$.
     */
    Scalar porosity() const
    { return porosity_; }

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


protected:

    Scalar porosity_; //!< Effective porosity within the control volume
    PermeabilityType permeability_; //!< Effective permeability within the control volume
    Scalar relativePermeability_[numPhases]; //!< Relative permeability within the control volume
    Scalar diffCoeff_[numPhases]; //!< Binary diffusion coefficients for the phases
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace Dumux

#endif
