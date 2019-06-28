// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>

#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>
#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mpadapter.hh>
#include "pressureformulation.hh"

namespace Dumux {

// forward declaration
template <class Traits, bool enableChemicalNonEquilibrium>
class MPNCVolumeVariablesImplementation;
/*!
 * \ingroup MPNCModel
 * \brief Contains the quantities which are constant within a finite volume in the MpNc model.
 *
 * \tparam Traits Class encapsulating types to be used by the vol vars
 */
template <class Traits>
using MPNCVolumeVariables =  MPNCVolumeVariablesImplementation<Traits, Traits::ModelTraits::enableChemicalNonEquilibrium()>;

template <class Traits>
class MPNCVolumeVariablesImplementation<Traits, false>
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, MPNCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, MPNCVolumeVariables<Traits> >;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;

    using ModelTraits = typename Traits::ModelTraits;
    static constexpr auto pressureFormulation = ModelTraits::pressureFormulation();

    static constexpr bool enableThermalNonEquilibrium = ModelTraits::enableThermalNonEquilibrium();
    static constexpr bool enableChemicalNonEquilibrium = ModelTraits::enableChemicalNonEquilibrium();
    static constexpr bool enableDiffusion = ModelTraits::enableMolecularDiffusion();

    using Indices = typename ModelTraits::Indices;
    using ComponentVector = Dune::FieldVector<Scalar, ModelTraits::numFluidComponents()>;
    using CompositionFromFugacities = Dumux::CompositionFromFugacities<Scalar, typename Traits::FluidSystem>;
    using EffDiffModel = typename Traits::EffectiveDiffusivityModel;

public:
    //! Export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export the fluid state type
    using FluidState = typename Traits::FluidState;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    //! Return number of phases considered by the model
    static constexpr int numFluidPhases() { return ModelTraits::numFluidPhases(); }
    //! Return number of components considered by the model
    static constexpr int numFluidComps = ParentType::numFluidComponents();

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        //calculate the remaining quantities
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        // relative permeabilities
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        using MPAdapter = MPAdapter<MaterialLaw, numFluidPhases()>;
        MPAdapter::relativePermeabilities(relativePermeability_, materialParams, fluidState_, wPhaseIdx);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        // porosity
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);

        if (enableDiffusion)
        {
            for (int phaseIdx = 0; phaseIdx < numFluidPhases(); ++phaseIdx)
            {
                int compIIdx = phaseIdx;
                for (unsigned int compJIdx = 0; compJIdx < numFluidComps; ++compJIdx)
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
                        setEffectiveDiffusionCoefficient_(phaseIdx, compJIdx);
                    }
                }
            }
        }

        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
        EnergyVolVars::updateEffectiveThermalConductivity();
    }

    /*!
     * \brief Sets complete fluid state.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param fluidState A container with the current (physical) state of the fluid
     * \param solidState A container with the current (physical) state of the solid
     */

    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {

        /////////////
        // set the fluid phase temperatures
        /////////////
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        /////////////
        // set the phase saturations
        /////////////
        auto&& priVars = elemSol[scv.localDofIndex()];
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numFluidPhases() - 1; ++phaseIdx) {
            sumSat += priVars[Indices::s0Idx + phaseIdx];
            fluidState.setSaturation(phaseIdx, priVars[Indices::s0Idx + phaseIdx]);
        }
        Valgrind::CheckDefined(sumSat);
        fluidState.setSaturation(numFluidPhases() - 1, 1.0 - sumSat);

        /////////////
        // set the phase pressures
        /////////////
        // capillary pressure parameters
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);
        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, scv, elemSol);
        // capillary pressures
        std::vector<Scalar> capPress(numFluidPhases());
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        using MPAdapter = MPAdapter<MaterialLaw, numFluidPhases()>;
        MPAdapter::capillaryPressures(capPress, materialParams, fluidState, wPhaseIdx);
        // add to the pressure of the first fluid phase

        // depending on which pressure is stored in the primary variables
        if(pressureFormulation == MpNcPressureFormulation::mostWettingFirst){
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            const Scalar pw = priVars[Indices::p0Idx];
            for (int phaseIdx = 0; phaseIdx < numFluidPhases(); ++phaseIdx)
                fluidState.setPressure(phaseIdx, pw - capPress[0] + capPress[phaseIdx]);
        }
        else if(pressureFormulation == MpNcPressureFormulation::leastWettingFirst){
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            const Scalar pn = priVars[Indices::p0Idx];
            for (int phaseIdx = numFluidPhases()-1; phaseIdx >= 0; --phaseIdx)
                fluidState.setPressure(phaseIdx, pn - capPress[numFluidPhases()-1] + capPress[phaseIdx]);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "MPNCVolumeVariables do not support the chosen pressure formulation");

        /////////////
        // set the fluid compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        ComponentVector fug;
        // retrieve component fugacities
        for (int compIdx = 0; compIdx < numFluidComps; ++compIdx)
            fug[compIdx] = priVars[Indices::fug0Idx + compIdx];

        // calculate phase compositions
        for (int phaseIdx = 0; phaseIdx < numFluidPhases(); ++phaseIdx) {
            // initial guess
            for (int compIdx = 0; compIdx < numFluidComps; ++compIdx) {
                Scalar x_ij = 1.0/numFluidComps;

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
        for (int phaseIdx = 0; phaseIdx < numFluidPhases(); ++phaseIdx) {
            // viscosities
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the enthalpy
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }

    /*!
     * \brief Returns the fluid configuration at the given primary
     *        variables.
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
     * \brief Returns the concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar molarity(const int phaseIdx, int compIdx) const
    { return fluidState_.molarity(phaseIdx, compIdx); }

    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the temperature inside the sub-control volume.
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
     * \brief Returns the internal energy \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
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
     * \brief Returns the fugacity \f$\mathrm{[kg/m^3]}\f$ the of the component.
     */
    Scalar fugacity(const int compIdx) const
    { return fluidState_.fugacity(compIdx); }

    /*!
     * \brief Returns the average molar mass \f$\mathrm{[kg/m^3]}\f$ the of the phase.
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
    { return solidState_.porosity(); }

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
     * \brief Returns the diffusion coefficient.
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
     * \brief Returns the diffusion coefficient.
     */
    Scalar effectiveDiffusivity(int phaseIdx, int compIdx) const
    {
        if (compIdx < phaseIdx)
            return effectiveDiffCoeff_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return effectiveDiffCoeff_[phaseIdx][compIdx-1];
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
        for (int compIdx = 0; compIdx < numFluidComps; ++compIdx)
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

    void setEffectiveDiffusionCoefficient_(int phaseIdx, int compIdx)
    {
        if (compIdx < phaseIdx)
            effectiveDiffCoeff_[phaseIdx][compIdx] = EffDiffModel::effectiveDiffusivity(*this, diffCoefficient_[phaseIdx][compIdx], phaseIdx);
        else if (compIdx > phaseIdx)
            effectiveDiffCoeff_[phaseIdx][compIdx-1] = EffDiffModel::effectiveDiffusivity(*this, diffCoefficient_[phaseIdx][compIdx-1], phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Effective diffusion coefficient for phaseIdx = compIdx doesn't exist");
    }

    std::array<std::array<Scalar, numFluidComps-1>, numFluidPhases()> diffCoefficient_;
    std::array<std::array<Scalar, numFluidComps-1>, numFluidPhases()> effectiveDiffCoeff_;
    Scalar porosity_; //!< Effective porosity within the control volume
    std::array<Scalar, ModelTraits::numFluidPhases()> relativePermeability_; //!< Effective relative permeability within the control volume
    PermeabilityType permeability_;

    //! Mass fractions of each component within each phase
    FluidState fluidState_;
    SolidState solidState_;

};

template <class Traits>
class MPNCVolumeVariablesImplementation<Traits, true>
    : public PorousMediumFlowVolumeVariables<Traits>
    , public EnergyVolumeVariables<Traits, MPNCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables< Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, MPNCVolumeVariables<Traits> >;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;

    using ModelTraits = typename Traits::ModelTraits;
    static constexpr auto pressureFormulation = ModelTraits::pressureFormulation();

    static constexpr bool enableThermalNonEquilibrium = ModelTraits::enableThermalNonEquilibrium();
    static constexpr bool enableChemicalNonEquilibrium = ModelTraits::enableChemicalNonEquilibrium();
    static constexpr bool enableDiffusion = ModelTraits::enableMolecularDiffusion();

    using Indices = typename ModelTraits::Indices;
    using ComponentVector = Dune::FieldVector<Scalar,  ModelTraits::numFluidComponents()>;
    using CompositionFromFugacities = Dumux::CompositionFromFugacities<Scalar, typename Traits::FluidSystem>;
    using EffDiffModel = typename Traits::EffectiveDiffusivityModel;

    using ParameterCache = typename Traits::FluidSystem::ParameterCache;
public:
    //! Export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export the fluid state type
    using FluidState = typename Traits::FluidState;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    //! Return number of phases considered by the model
    static constexpr int numFluidPhases() { return ModelTraits::numFluidPhases(); }
    //! Return number of components considered by the model
    static constexpr int numFluidComps = ParentType::numFluidComponents();
    using ConstraintSolver = MiscibleMultiPhaseComposition<Scalar, FluidSystem>;

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        // calculate the remaining quantities
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);

        // relative permeabilities
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        using MPAdapter = MPAdapter<MaterialLaw, numFluidPhases()>;
        MPAdapter::relativePermeabilities(relativePermeability_,  materialParams, fluidState_, wPhaseIdx);
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);

        if (enableDiffusion)
        {
            for (int phaseIdx = 0; phaseIdx < numFluidPhases(); ++phaseIdx)
            {
                int compIIdx = phaseIdx;
                for (unsigned int compJIdx = 0; compJIdx < numFluidComps; ++compJIdx)
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
                        setEffectiveDiffusionCoefficient_(phaseIdx, compJIdx);
                    }
                }
            }
        }

        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
        EnergyVolVars::updateEffectiveThermalConductivity();
    }

    /*!
     * \brief Sets complete fluid state.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param fluidState A container with the current (physical) state of the fluid
     * \param solidState A container with the current (physical) state of the solid
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {
        /////////////
        // set the fluid phase temperatures
        /////////////
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);
        /////////////
        // set the phase saturations
        /////////////
        auto&& priVars = elemSol[scv.localDofIndex()];
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numFluidPhases() - 1; ++phaseIdx) {
            sumSat += priVars[Indices::s0Idx + phaseIdx];
            fluidState.setSaturation(phaseIdx, priVars[Indices::s0Idx + phaseIdx]);
        }
        Valgrind::CheckDefined(sumSat);
        fluidState.setSaturation(numFluidPhases() - 1, 1.0 - sumSat);

        /////////////
        // set the phase pressures
        /////////////
        // capillary pressure parameters
        const auto& materialParams =
            problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);
        // capillary pressures
        std::vector<Scalar> capPress(numFluidPhases());
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        using MPAdapter = MPAdapter<MaterialLaw, numFluidPhases()>;
        MPAdapter::capillaryPressures(capPress, materialParams, fluidState, wPhaseIdx);
        // add to the pressure of the first fluid phase

        // depending on which pressure is stored in the primary variables
        if(pressureFormulation == MpNcPressureFormulation::mostWettingFirst){
            // This means that the pressures are sorted from the most wetting to the least wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pw
            const Scalar pw = priVars[Indices::p0Idx];
            for (int phaseIdx = 0; phaseIdx < numFluidPhases(); ++phaseIdx)
                fluidState.setPressure(phaseIdx, pw - capPress[0] + capPress[phaseIdx]);
        }
        else if(pressureFormulation == MpNcPressureFormulation::leastWettingFirst){
            // This means that the pressures are sorted from the least wetting to the most wetting-1 in the primary variables vector.
            // For two phases this means that there is one pressure as primary variable: pn
            const Scalar pn = priVars[Indices::p0Idx];
            for (int phaseIdx = numFluidPhases()-1; phaseIdx >= 0; --phaseIdx)
                fluidState.setPressure(phaseIdx, pn - capPress[numFluidPhases()-1] + capPress[phaseIdx]);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "MPNCVolumeVariables do not support the chosen pressure formulation");


        /////////////
        // set the fluid compositions
        /////////////
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        ComponentVector fug;
        // retrieve component fugacities
        for (int compIdx = 0; compIdx < numFluidComps; ++compIdx)
            fug[compIdx] = priVars[Indices::fug0Idx + compIdx];

         updateMoleFraction(fluidState,
                            paramCache,
                            priVars);


        // dynamic viscosities
        for (int phaseIdx = 0; phaseIdx < numFluidPhases(); ++phaseIdx) {
            // viscosities
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            fluidState.setViscosity(phaseIdx, mu);

            // compute and set the enthalpy
            Scalar h = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
            fluidState.setEnthalpy(phaseIdx, h);
        }
   }

    /*!
     * \brief Updates composition of all phases in the mutable
     *        parameters from the primary variables.
     *
     *        \param actualFluidState Container for all the secondary variables concerning the fluids
     *        \param paramCache Container for cache parameters
     *        \param priVars The primary Variables
     */
    void updateMoleFraction(FluidState & actualFluidState,
                            ParameterCache & paramCache,
                            const typename Traits::PrimaryVariables& priVars)
    {
        // setting the mole fractions of the fluid state
        for(int phaseIdx=0; phaseIdx<numFluidPhases(); ++phaseIdx)
        {
                // set the component mole fractions
                for (int compIdx = 0; compIdx < numFluidComps; ++compIdx) {
                    actualFluidState.setMoleFraction(phaseIdx,
                           compIdx,
                           priVars[Indices::moleFrac00Idx +
                                   phaseIdx*numFluidComps +
                                   compIdx]);
                }
        }

//          // For using the ... other way of calculating equilibrium
//          THIS IS ONLY FOR silencing Valgrind but is not used in this model
        for(int phaseIdx=0; phaseIdx<numFluidPhases(); ++phaseIdx)
            for (int compIdx = 0; compIdx < numFluidComps; ++compIdx) {
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
                                    paramCache) ;

        // Setting the equilibrium composition (in a kinetic model not necessarily the same as the actual mole fraction)
        for(int phaseIdx=0; phaseIdx<numFluidPhases(); ++phaseIdx){
            for (int compIdx=0; compIdx< numFluidComps; ++ compIdx){
                xEquil_[phaseIdx][compIdx] = equilFluidState.moleFraction(phaseIdx, compIdx);
            }
        }

        // compute densities of all phases
        for(int phaseIdx=0; phaseIdx<numFluidPhases(); ++phaseIdx){
            const Scalar rho = FluidSystem::density(actualFluidState, paramCache, phaseIdx);
            actualFluidState.setDensity(phaseIdx, rho);
            const Scalar rhoMolar = FluidSystem::molarDensity(actualFluidState, paramCache, phaseIdx);
            actualFluidState.setMolarDensity(phaseIdx, rhoMolar);
        }

    }

    /*!
     * \brief The mole fraction we would have in the case of chemical equilibrium /
     *        on the interface.
     *
     * \param phaseIdx The index of the fluid phase
     * \param compIdx The local index of the component
     */
    const Scalar xEquil(const unsigned int phaseIdx, const unsigned int compIdx) const
    {
        return xEquil_[phaseIdx][compIdx] ;
    }

    /*!
     * \brief Returns the fluid configuration at the given primary
     *        variables
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
     * \brief Returns the concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     */
    Scalar molarity(const int phaseIdx, int compIdx) const
    { return fluidState_.molarity(phaseIdx, compIdx); }

    /*!
     * \brief Returns the molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    { return fluidState_.molarDensity(phaseIdx);}

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(const int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(const int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(0/* phaseIdx*/); }

    /*!
     * \brief Returns the enthalpy \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     */
    Scalar enthalpy(const int phaseIdx) const
    { return fluidState_.enthalpy(phaseIdx); }

    /*!
     * \brief Returns the internal energy \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
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
     * \brief Returns the fugacity \f$\mathrm{[kg/m^3]}\f$ the of the component.
     */
    Scalar fugacity(const int compIdx) const
    { return fluidState_.fugacity(compIdx); }

    /*!
     * \brief Returns the average molar mass \f$\mathrm{[kg/m^3]}\f$ the of the phase.
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
    { return solidState_.porosity(); }

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
     * \brief Returns the effective diffusion coefficients for a phase in \f$[m^2/s]\f$.
     */
    Scalar effectiveDiffusivity(int phaseIdx, int compIdx) const
    {
        if (compIdx < phaseIdx)
            return effectiveDiffCoeff_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            return effectiveDiffCoeff_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient called for phaseIdx = compIdx");
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
        for (int compIdx = 0; compIdx < numFluidComps; ++compIdx)
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

    void setEffectiveDiffusionCoefficient_(int phaseIdx, int compIdx)
    {
        if (compIdx < phaseIdx)
            effectiveDiffCoeff_[phaseIdx][compIdx] = EffDiffModel::effectiveDiffusivity(*this, diffCoefficient_[phaseIdx][compIdx], phaseIdx);
        else if (compIdx > phaseIdx)
            effectiveDiffCoeff_[phaseIdx][compIdx-1] = EffDiffModel::effectiveDiffusivity(*this, diffCoefficient_[phaseIdx][compIdx-1], phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Effective diffusion coefficient for phaseIdx = compIdx doesn't exist");
    }

    std::array<std::array<Scalar, numFluidComps-1>, numFluidPhases()> diffCoefficient_;
    std::array<std::array<Scalar, numFluidComps-1>, numFluidPhases()> effectiveDiffCoeff_;
    std::array<Scalar, ModelTraits::numFluidPhases()> relativePermeability_; //!< Effective relative permeability within the control volume
    PermeabilityType permeability_;
    std::array<std::array<Scalar, numFluidComps>, numFluidPhases()> xEquil_;

    //! Mass fractions of each component within each phase
    FluidState fluidState_;
    SolidState solidState_;
};

} // end namespace

#endif
