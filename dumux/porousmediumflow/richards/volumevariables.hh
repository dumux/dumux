// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsModel
 * \brief Volume averaged quantities required by the Richards model.
 */

#ifndef DUMUX_RICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_VOLUME_VARIABLES_HH

#include <cassert>

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/state.hh>

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/material/idealgas.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Volume averaged quantities required by the Richards model.
 *
 * This contains the quantities which are are constant within a finite
 * volume in the Richards model.
 */
template <class Traits>
class RichardsVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, RichardsVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, RichardsVolumeVariables<Traits> >;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using ModelTraits = typename Traits::ModelTraits;

    static constexpr int numFluidComps = ParentType::numFluidComponents();
    static constexpr int numPhases = ParentType::numFluidPhases();

    // checks if the fluid system uses the Richards model index convention
    static constexpr auto fsCheck = ModelTraits::checkFluidSystem(typename Traits::FluidSystem{});

public:
    //! Export type of the fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of the fluid state
    using FluidState = typename Traits::FluidState;
    //! Export type of the fluid state
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    using Indices = typename Traits::ModelTraits::Indices;

    //! Export phase indices
    static constexpr auto liquidPhaseIdx = Traits::FluidSystem::phase0Idx;
    static constexpr auto gasPhaseIdx = Traits::FluidSystem::phase1Idx;

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
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);

        // precompute the minimum capillary pressure (entry pressure)
        // needed to make sure we don't compute unphysical capillary pressures and thus saturations
        minPc_ = fluidMatrixInteraction.endPointPc();

        //update porosity before calculating the effective properties depending on it
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        //////////
        // specify the other parameters
        //////////
        relativePermeabilityWetting_ = fluidMatrixInteraction.krw(fluidState_.saturation(liquidPhaseIdx));
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
        EnergyVolVars::updateEffectiveThermalConductivity();
    }

    /*!
     * \brief Fills the fluid state according to the primary variables.
     *
     * Taking the information from the primary variables,
     * the fluid state is filled with every information that is
     * necessary to evaluate the model's local residual.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The problem at hand.
     * \param element The current element.
     * \param scv The subcontrol volume.
     * \param fluidState The fluid state to fill.
     * \param solidState The solid state to fill.
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

        const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);

        const auto& priVars = elemSol[scv.localDofIndex()];

        // set the wetting pressure
        using std::max;
        fluidState.setPressure(liquidPhaseIdx, priVars[Indices::pressureIdx]);
        fluidState.setPressure(gasPhaseIdx, max(problem.nonwettingReferencePressure(), fluidState.pressure(liquidPhaseIdx) + minPc_));

        // compute the capillary pressure to compute the saturation
        // make sure that we the capillary pressure is not smaller than the minimum pc
        // this would possibly return unphysical values from regularized material laws
        using std::max;
        const Scalar pc = max(minPc_, problem.nonwettingReferencePressure() - fluidState.pressure(liquidPhaseIdx));
        const Scalar sw = fluidMatrixInteraction.sw(pc);
        fluidState.setSaturation(liquidPhaseIdx, sw);
        fluidState.setSaturation(gasPhaseIdx, 1.0-sw);

        // density and viscosity
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        fluidState.setDensity(liquidPhaseIdx,
                              FluidSystem::density(fluidState, paramCache, liquidPhaseIdx));
        fluidState.setDensity(gasPhaseIdx,
                              FluidSystem::density(fluidState, paramCache, gasPhaseIdx));

        fluidState.setViscosity(liquidPhaseIdx,
                                FluidSystem::viscosity(fluidState, paramCache, liquidPhaseIdx));

        // compute and set the enthalpy
        fluidState.setEnthalpy(liquidPhaseIdx, EnergyVolVars::enthalpy(fluidState, paramCache, liquidPhaseIdx));
        fluidState.setEnthalpy(gasPhaseIdx, EnergyVolVars::enthalpy(fluidState, paramCache, gasPhaseIdx));
    }

    /*!
     * \brief Returns the fluid configuration at the given primary
     *        variables.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the temperature.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the average porosity [] within the control volume.
     *
     * The porosity is defined as the ratio of the pore space to the
     * total volume, i.e. \f[ \Phi := \frac{V_{pore}}{V_{pore} + V_{rock}} \f]
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

    /*!
     * \brief Returns the average absolute saturation [] of a given
     *        fluid phase within the finite volume.
     *
     * The saturation of a fluid phase is defined as the fraction of
     * the pore volume filled by it, i.e.
     * \f[ S_\alpha := \frac{V_\alpha}{V_{pore}} = \phi \frac{V_\alpha}{V} \f]
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar saturation(const int phaseIdx = liquidPhaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the average mass density \f$\mathrm{[kg/m^3]}\f$ of a given
     *        fluid phase within the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar density(const int phaseIdx = liquidPhaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * For the nonwetting phase (i.e. the gas phase), we assume
     * infinite mobility, which implies that the nonwetting phase
     * pressure is equal to the finite volume's reference pressure
     * defined by the problem.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar pressure(const int phaseIdx = liquidPhaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns the effective mobility \f$\mathrm{[1/(Pa*s)]}\f$ of a given phase within
     *        the control volume.
     *
     * The mobility of a fluid phase is defined as the relative
     * permeability of the phase (given by the chosen material law)
     * divided by the dynamic viscosity of the fluid, i.e.
     * \f[ \lambda_\alpha := \frac{k_{r\alpha}}{\mu_\alpha} \f]
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar mobility(const int phaseIdx = liquidPhaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     * \note The nonwetting phase is infinitely mobile
     */
    Scalar viscosity(const int phaseIdx = liquidPhaseIdx) const
    { return phaseIdx == liquidPhaseIdx ? fluidState_.viscosity(liquidPhaseIdx) : 0.0; }

    /*!
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar relativePermeability(const int phaseIdx = liquidPhaseIdx) const
    { return phaseIdx == liquidPhaseIdx ? relativePermeabilityWetting_ : 1.0; }

    /*!
     * \brief Returns the effective capillary pressure \f$\mathrm{[Pa]}\f$ within the
     *        control volume.
     *
     * The capillary pressure is defined as the difference in
     * pressures of the nonwetting and the wetting phase, i.e.
     * \f[ p_c = p_n - p_w \f]
     *
     * \note Capillary pressures are always larger than the entry pressure
     *       This regularization doesn't affect the residual in which pc is not needed.
     */
    Scalar capillaryPressure() const
    {
        using std::max;
        return max(minPc_, pressure(gasPhaseIdx) - pressure(liquidPhaseIdx));
    }

    /*!
     * \brief Returns the pressureHead \f$\mathrm{[cm]}\f$ of a given phase within
     *        the control volume.
     *
     * For the nonwetting phase (i.e. the gas phase), we assume
     * infinite mobility, which implies that the nonwetting phase
     * pressure is equal to the finite volume's reference pressure
     * defined by the problem.
     *
     * \param phaseIdx The index of the fluid phase
     * \note this function is here as a convenience to the user to not have to
     *       manually do a conversion. It is not correct if the density is not constant
     *       or the gravity different
     */
    Scalar pressureHead(const int phaseIdx = liquidPhaseIdx) const
    { return 100.0 *(pressure(phaseIdx) - pressure(gasPhaseIdx))/density(phaseIdx)/9.81; }

    /*!
     * \brief Returns the water content of a fluid phase within the finite volume.
     *
     * The water content is defined as the fraction of
     * the saturation divided by the porosity.

     * \param phaseIdx The index of the fluid phase
     * \note this function is here as a convenience to the user to not have to
     *       manually do a conversion.
     */
    Scalar waterContent(const int phaseIdx = liquidPhaseIdx) const
    { return saturation(phaseIdx) * solidState_.porosity(); }

protected:
    FluidState fluidState_; //!< the fluid state
    SolidState solidState_;
    Scalar relativePermeabilityWetting_; //!< the relative permeability of the wetting phase
    PermeabilityType permeability_; //!< the intrinsic permeability
    Scalar minPc_; //!< the minimum capillary pressure (entry pressure)

    // Effective diffusion coefficients for the phases
    Scalar effectiveDiffCoeff_;

};
} // end namespace Dumux

#endif
