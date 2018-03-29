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
 * \ingroup MPNCModel
 * \brief Contains the secondary variables (Quantities which are
 *        constant within a finite volume) of the MpNc model.
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_HH

#include "indices.hh"

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/volumevariables.hh>

#include <dumux/material/constraintsolvers/ncpflash.hh>
#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>

namespace Dumux
{

// forward declaration
template <class TypeTag, bool enableChemicalNonEquilibrium>
class MPNCVolumeVariablesImplementation;
/*!
 * \ingroup MPNCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the MpNc model.
 */

template <class TypeTag>
using MPNCVolumeVariables = MPNCVolumeVariablesImplementation<TypeTag, GET_PROP_TYPE(TypeTag, ModelTraits)::enableChemicalNonEquilibrium()>;

template <class TypeTag>
class MPNCVolumeVariablesImplementation<TypeTag, false>
    : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using Element = typename GridView::template Codim<0>::Entity;

    // formulations
    enum {
        pressureFormulation = GET_PROP_VALUE(TypeTag, PressureFormulation),
        mostWettingFirst    = MpNcPressureFormulation::mostWettingFirst,
        leastWettingFirst   = MpNcPressureFormulation::leastWettingFirst
    };

    enum {numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases()};
    enum {numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents()};
    enum {s0Idx = Indices::s0Idx};
    enum {p0Idx = Indices::p0Idx};
    enum {fug0Idx = Indices::fug0Idx};

    static constexpr bool enableDiffusion = GET_PROP_TYPE(TypeTag, ModelTraits)::enableMolecularDiffusion();


    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using CompositionFromFugacities = Dumux::CompositionFromFugacities<Scalar, FluidSystem>;

public:
    // export type of fluid state for non-isothermal models
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParentType::ParentType;

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_);

        //calculate the remaining quantities
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        // relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams,
                                            fluidState_);
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);
        if (enableDiffusion)
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                int compIIdx = phaseIdx;
                for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
                {
                    // binary diffusion coefficients
                    if(compIIdx!= compJIdx)
                    {
                        setDiffusionCoefficient_(phaseIdx, compJIdx,
                                                FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                        paramCache,
                                                                                        phaseIdx,
                                                                                        compIIdx,
                                                                                        compJIdx));
                    }
                }
            }
        }

        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     * \param isOldSol Specifies whether this is the previous solution or the current one
     */
     template<class ElementSolution>
     void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        /////////////
        // set the fluid phase temperatures
        /////////////
        ParentType::updateTemperature(elemSol, problem, element, scv, fluidState);

        /////////////
        // set the phase saturations
        /////////////
        auto&& priVars = ParentType::extractDofPriVars(elemSol, scv);
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            sumSat += priVars[s0Idx + phaseIdx];
            fluidState.setSaturation(phaseIdx, priVars[s0Idx + phaseIdx]);
        }
        Valgrind::CheckDefined(sumSat);
        fluidState.setSaturation(numPhases - 1, 1.0 - sumSat);

        /////////////
        // set the phase pressures
        /////////////
        // capillary pressure parameters
        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, scv, elemSol);
        // capillary pressures
        Scalar capPress[numPhases];
        MaterialLaw::capillaryPressures(capPress, materialParams, fluidState);
        // add to the pressure of the first fluid phase

        // depending on which pressure is stored in the primary variables
        if(pressureFormulation == mostWettingFirst){
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            const Scalar pw = priVars[p0Idx];
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx, pw - capPress[0] + capPress[phaseIdx]);
        }
        else if(pressureFormulation == leastWettingFirst){
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            const Scalar pn = priVars[p0Idx];
            for (int phaseIdx = numPhases-1; phaseIdx >= 0; --phaseIdx)
                fluidState.setPressure(phaseIdx, pn - capPress[numPhases-1] + capPress[phaseIdx]);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << pressureFormulation << " is invalid.");

        /////////////
        // set the fluid compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        ComponentVector fug;
        // retrieve component fugacities
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            fug[compIdx] = priVars[fug0Idx + compIdx];

            // calculate phase compositions
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // initial guess
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar x_ij = 1.0/numComponents;

                    // set initial guess of the component's mole fraction
                    fluidState.setMoleFraction(phaseIdx,
                                            compIdx,
                                            x_ij);
                }
                // calculate the phase composition from the component
                // fugacities
                CompositionFromFugacities::guessInitial(fluidState, paramCache, phaseIdx, fug);
                CompositionFromFugacities::solve(fluidState, paramCache, phaseIdx, fug);
            }
        // dynamic viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // viscosities
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the enthalpy
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }

        // make sure the quantities in the fluid state are well-defined
        fluidState.checkDefined();
   }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
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
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar molarity(const int phaseIdx, int compIdx) const
    { return fluidState_.molarity(phaseIdx, compIdx); }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(0/* phaseIdx*/); }

    Scalar temperature(const int phaseIdx) const
    { return fluidState_.temperature(phaseIdx); }

    /*!
     * \brief Return enthalpy \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     */
    Scalar enthalpy(const int phaseIdx) const
    { return fluidState_.enthalpy(phaseIdx); }

    /*!
     * \brief Return internal energy \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     */
    Scalar internalEnergy(const int phaseIdx) const
    { return fluidState_.internalEnergy(phaseIdx); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of a fluid phase in
     *        the sub-control volume.
     */
    Scalar fluidThermalConductivity(const int phaseIdx) const
    { return FluidSystem::thermalConductivity(fluidState_, phaseIdx); }

    /*!
     * \brief Return fugacity \f$\mathrm{[kg/m^3]}\f$ the of the component.
     */
    Scalar fugacity(const int compIdx) const
    { return fluidState_.fugacity(compIdx); }

    /*!
     * \brief Return average molar mass \f$\mathrm{[kg/m^3]}\f$ the of the phase.
     */
    Scalar averageMolarMass(const int phaseIdx) const
    { return fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     *        \param phaseIdx The local index of the phases
     */
    Scalar mobility(const unsigned int phaseIdx) const
    {
        return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Returns the viscosity of a given phase within
     *        the control volume.
     */
    Scalar viscosity(const unsigned int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the relative permeability of a given phase within
     *        the control volume.
     *
     *       \param phaseIdx The local index of the phases
     */
    Scalar relativePermeability(const unsigned int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns true if the fluid state is in the active set
     *        for a phase,
     *
     *        \param phaseIdx The local index of the phases
     */
    bool isPhaseActive(const unsigned int phaseIdx) const
    {
        return
            phasePresentIneq(fluidState(), phaseIdx) -
            phaseNotPresentIneq(fluidState(), phaseIdx)
            >= 0;
    }

    /*!
     * \brief Returns the diffusion coefficient
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        if (compIdx < phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient called for phaseIdx = compIdx");
    }

    /*!
     * \brief Returns the value of the NCP-function for a phase.
     *
     *      \param phaseIdx The local index of the phases
     */
    Scalar phaseNcp(const unsigned int phaseIdx) const
    {
        Scalar aEval = phaseNotPresentIneq(fluidState(), phaseIdx);
        Scalar bEval = phasePresentIneq(fluidState(), phaseIdx);
        if (aEval > bEval)
            return phasePresentIneq(fluidState(), phaseIdx);
        return phaseNotPresentIneq(fluidState(), phaseIdx);
    };

    /*!
     * \brief Returns the value of the inequality where a phase is
     *        present.
     *
     *        \param phaseIdx The local index of the phases
     *        \param fluidState Container for all the secondary variables concerning the fluids
     */
    Scalar phasePresentIneq(const FluidState &fluidState,
                            const unsigned int phaseIdx) const
    { return fluidState.saturation(phaseIdx); }

    /*!
     * \brief Returns the value of the inequality where a phase is not
     *        present.
     *
     *        \param phaseIdx The local index of the phases
     *        \param fluidState Container for all the secondary variables concerning the fluids
     */
    Scalar phaseNotPresentIneq(const FluidState &fluidState,
                               const unsigned int phaseIdx) const
    {
        // difference of sum of mole fractions in the phase from 100%
        Scalar a = 1;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            a -= fluidState.moleFraction(phaseIdx, compIdx);
        return a;
    }

protected:

    void setDiffusionCoefficient_(int phaseIdx, int compIdx, Scalar d)
    {
        if (compIdx < phaseIdx)
            diffCoefficient_[phaseIdx][compIdx] = std::move(d);
        else if (compIdx > phaseIdx)
            diffCoefficient_[phaseIdx][compIdx-1] = std::move(d);
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient for phaseIdx = compIdx doesn't exist");
    }

    std::array<std::array<Scalar, numComponents-1>, numPhases> diffCoefficient_;
    Scalar porosity_; //!< Effective porosity within the control volume
    Scalar relativePermeability_[numPhases]; //!< Effective relative permeability within the control volume
    PermeabilityType permeability_;

    //! Mass fractions of each component within each phase
    FluidState fluidState_;
};

template <class TypeTag>
class MPNCVolumeVariablesImplementation<TypeTag, true>
     : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParameterCache = typename FluidSystem::ParameterCache;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);

    // formulations
    enum {
        pressureFormulation = GET_PROP_VALUE(TypeTag, PressureFormulation),
        mostWettingFirst    = MpNcPressureFormulation::mostWettingFirst,
        leastWettingFirst   = MpNcPressureFormulation::leastWettingFirst
    };

    enum {numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases()};
    enum {numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents()};
    enum {s0Idx = Indices::s0Idx};
    enum {p0Idx = Indices::p0Idx};
    enum {fug0Idx = Indices::fug0Idx};
    enum { moleFrac00Idx = Indices::conti0EqIdx };

    static constexpr bool enableDiffusion = GET_PROP_TYPE(TypeTag, ModelTraits)::enableMolecularDiffusion();

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using ConstraintSolver = MiscibleMultiPhaseComposition<Scalar, FluidSystem>;

public:
    // export type of fluid state for non-isothermal models
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using ParentType::ParentType;

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_);

        //calculate the remaining quantities
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        // relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams,
                                            fluidState_);
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);
        if (enableDiffusion)
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                int compIIdx = phaseIdx;
                for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
                {
                    // binary diffusion coefficients
                    if(compIIdx!= compJIdx)
                    {
                        setDiffusionCoefficient_(phaseIdx, compJIdx,
                                                FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                        paramCache,
                                                                                        phaseIdx,
                                                                                        compIIdx,
                                                                                        compJIdx));
                    }
                }
            }
        }

        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     * \param isOldSol Specifies whether this is the previous solution or the current one
     */
     template<class ElementSolution>
     void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        /////////////
        // set the fluid phase temperatures
        /////////////
        ParentType::updateTemperature(elemSol, problem, element, scv, fluidState);
        /////////////
        // set the phase saturations
        /////////////
        auto&& priVars = ParentType::extractDofPriVars(elemSol, scv);
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            sumSat += priVars[s0Idx + phaseIdx];
            fluidState.setSaturation(phaseIdx, priVars[s0Idx + phaseIdx]);
        }
        Valgrind::CheckDefined(sumSat);
        fluidState.setSaturation(numPhases - 1, 1.0 - sumSat);

        /////////////
        // set the phase pressures
        /////////////
        // capillary pressure parameters
        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, scv, elemSol);
        // capillary pressures
        Scalar capPress[numPhases];
        MaterialLaw::capillaryPressures(capPress, materialParams, fluidState);
        // add to the pressure of the first fluid phase

        // depending on which pressure is stored in the primary variables
        if(pressureFormulation == mostWettingFirst){
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            const Scalar pw = priVars[p0Idx];
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx, pw - capPress[0] + capPress[phaseIdx]);
        }
        else if(pressureFormulation == leastWettingFirst){
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            const Scalar pn = priVars[p0Idx];
            for (int phaseIdx = numPhases-1; phaseIdx >= 0; --phaseIdx)
                fluidState.setPressure(phaseIdx, pn - capPress[numPhases-1] + capPress[phaseIdx]);
        }
        else DUNE_THROW(Dune::InvalidStateException, "Formulation: " << pressureFormulation << " is invalid.");

        /////////////
        // set the fluid compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        ComponentVector fug;
        // retrieve component fugacities
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            fug[compIdx] = priVars[fug0Idx + compIdx];

         updateMoleFraction(fluidState,
                            paramCache,
                            priVars);


        // dynamic viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // viscosities
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the enthalpy
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }

        // make sure the quantities in the fluid state are well-defined
        fluidState.checkDefined();
   }

        /*!
     * \brief Update composition of all phases in the mutable
     *        parameters from the primary variables.
     *
     *        \param actualFluidState Container for all the secondary variables concerning the fluids
     *        \param paramCache Container for cache parameters
     *        \param priVars The primary Variables
     *        \param *hint the volume variables, usable for initial guess of composition
     */
    void updateMoleFraction(FluidState & actualFluidState,
                            ParameterCache & paramCache,
                            const PrimaryVariables & priVars)
    {
        // setting the mole fractions of the fluid state
        for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
        {
                // set the component mole fractions
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    actualFluidState.setMoleFraction(phaseIdx,
                           compIdx,
                           priVars[moleFrac00Idx +
                                   phaseIdx*numComponents +
                                   compIdx]);
                }
        }

//          // For using the ... other way of calculating equilibrium
//          THIS IS ONLY FOR silencing Valgrind but is not used in this model
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Scalar phi = FluidSystem::fugacityCoefficient(actualFluidState,
                                                                        paramCache,
                                                                        phaseIdx,
                                                                        compIdx);
                    actualFluidState.setFugacityCoefficient(phaseIdx,
                                                      compIdx,
                                                      phi);
            }

            FluidState equilFluidState; // the fluidState *on the interface* i.e. chemical equilibrium
            equilFluidState.assign(actualFluidState) ;
            ConstraintSolver::solve(equilFluidState,
                                    paramCache,
                                    /*setViscosity=*/false,
                                    /*setEnthalpy=*/false) ;

            // Setting the equilibrium composition (in a kinetic model not necessarily the same as the actual mole fraction)
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                for (int compIdx=0; compIdx< numComponents; ++ compIdx){
                    xEquil_[phaseIdx][compIdx] = equilFluidState.moleFraction(phaseIdx, compIdx);
                }
            }

            // compute densities of all phases
            for(int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                const Scalar rho = FluidSystem::density(actualFluidState, paramCache, phaseIdx);
                actualFluidState.setDensity(phaseIdx, rho);
            }

        }

    /*!
     * \brief The mole fraction we would have in the case of chemical equilibrium /
     *        on the interface.
     *
     *     \param phaseIdx The index of the fluid phase
     *     \param compIdx The local index of the component
     */
    const Scalar xEquil(const unsigned int phaseIdx, const unsigned int compIdx) const
    {
        return xEquil_[phaseIdx][compIdx] ;
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
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
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar molarity(const int phaseIdx, int compIdx) const
    { return fluidState_.molarity(phaseIdx, compIdx); }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    {
        return fluidState_.temperature(0/* phaseIdx*/);
    }

    Scalar temperature(const int phaseIdx) const
    { return fluidState_.temperature(phaseIdx); }

    /*!
     * \brief Return enthalpy \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     */
    Scalar enthalpy(const int phaseIdx) const
    { return fluidState_.enthalpy(phaseIdx); }

    /*!
     * \brief Return internal energy \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     */
    Scalar internalEnergy(const int phaseIdx) const
    { return fluidState_.internalEnergy(phaseIdx); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of a fluid phase in
     *        the sub-control volume.
     */
    Scalar fluidThermalConductivity(const int phaseIdx) const
    { return FluidSystem::thermalConductivity(fluidState_, phaseIdx); }

    /*!
     * \brief Return fugacity \f$\mathrm{[kg/m^3]}\f$ the of the component.
     */
    Scalar fugacity(const int compIdx) const
    { return fluidState_.fugacity(compIdx); }

    /*!
     * \brief Return average molar mass \f$\mathrm{[kg/m^3]}\f$ the of the phase.
     */
    Scalar averageMolarMass(const int phaseIdx) const
    { return fluidState_.averageMolarMass(phaseIdx); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     *        \param phaseIdx The local index of the phases
     */
    Scalar mobility(const unsigned int phaseIdx) const
    {
        return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Returns the viscosity of a given phase within
     *        the control volume.
     */
    Scalar viscosity(const unsigned int phaseIdx) const
    { return fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the relative permeability of a given phase within
     *        the control volume.
     *
     *       \param phaseIdx The local index of the phases
     */
    Scalar relativePermeability(const unsigned int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns true if the fluid state is in the active set
     *        for a phase,
     *
     *        \param phaseIdx The local index of the phases
     */
    bool isPhaseActive(const unsigned int phaseIdx) const
    {
        return
            phasePresentIneq(fluidState(), phaseIdx) -
            phaseNotPresentIneq(fluidState(), phaseIdx)
            >= 0;
    }

    /*!
     * \brief Returns the diffusion coefficient
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        if (compIdx < phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return diffCoefficient_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient called for phaseIdx = compIdx");
    }

    /*!
     * \brief Returns the value of the NCP-function for a phase.
     *
     *      \param phaseIdx The local index of the phases
     */
    Scalar phaseNcp(const unsigned int phaseIdx) const
    {
        Scalar aEval = phaseNotPresentIneq(fluidState(), phaseIdx);
        Scalar bEval = phasePresentIneq(fluidState(), phaseIdx);
        if (aEval > bEval)
            return phasePresentIneq(fluidState(), phaseIdx);
        return phaseNotPresentIneq(fluidState(), phaseIdx);
    };

    /*!
     * \brief Returns the value of the inequality where a phase is
     *        present.
     *
     *        \param phaseIdx The local index of the phases
     *        \param fluidState Container for all the secondary variables concerning the fluids
     */
    Scalar phasePresentIneq(const FluidState &fluidState,
                            const unsigned int phaseIdx) const
    { return fluidState.saturation(phaseIdx); }

    /*!
     * \brief Returns the value of the inequality where a phase is not
     *        present.
     *
     *        \param phaseIdx The local index of the phases
     *        \param fluidState Container for all the secondary variables concerning the fluids
     */
    Scalar phaseNotPresentIneq(const FluidState &fluidState,
                               const unsigned int phaseIdx) const
    {
        // difference of sum of mole fractions in the phase from 100%
        Scalar a = 1;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            a -= fluidState.moleFraction(phaseIdx, compIdx);
        return a;
    }

protected:

    void setDiffusionCoefficient_(int phaseIdx, int compIdx, Scalar d)
    {
        if (compIdx < phaseIdx)
            diffCoefficient_[phaseIdx][compIdx] = std::move(d);
        else if (compIdx > phaseIdx)
            diffCoefficient_[phaseIdx][compIdx-1] = std::move(d);
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient for phaseIdx = compIdx doesn't exist");
    }

    std::array<std::array<Scalar, numComponents-1>, numPhases> diffCoefficient_;
    Scalar porosity_; //!< Effective porosity within the control volume
    Scalar relativePermeability_[numPhases]; //!< Effective relative permeability within the control volume
    PermeabilityType permeability_;
    Scalar xEquil_[numPhases][numComponents];

    //! Mass fractions of each component within each phase
    FluidState fluidState_;
};

} // end namespace

#endif
