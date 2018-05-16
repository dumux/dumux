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
 * \ingroup RichardsNCModel
 * \brief  Contains the quantities which are constant within a
 *        finite volume in the Richards, n-component model.
 */
#ifndef DUMUX_RICHARDSNC_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDSNC_VOLUME_VARIABLES_HH

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

namespace Dumux {

/*!
 * \ingroup RichardsNCModel
 * \brief  Contains the quantities which are constant within a
 *        finite volume in the Richards, n-component model.
 */
template <class Traits>
class RichardsNCVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, RichardsNCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, RichardsNCVolumeVariables<Traits> >;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using Idx = typename Traits::ModelTraits::Indices;

    static constexpr int fluidSystemPhaseIdx = Idx::fluidSystemPhaseIdx;
    static constexpr bool useMoles = Traits::ModelTraits::useMoles();

public:
    //! export type of the fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! export type of the fluid state
    using FluidState = typename Traits::FluidState;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! export indices
    using Indices = typename Traits::ModelTraits::Indices;
    //! export phase acess indices
    static constexpr int liquidPhaseIdx = fluidSystemPhaseIdx;
    static constexpr int gasPhaseIdx = 1 - liquidPhaseIdx;

    /*!
     * \brief Update all quantities for a given control volume
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

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);
        //////////
        // specify the other parameters
        //////////
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        relativePermeabilityWetting_ = MaterialLaw::krw(materialParams, fluidState_.saturation(fluidSystemPhaseIdx));

        // precompute the minimum capillary pressure (entry pressure)
        // needed to make sure we don't compute unphysical capillary pressures and thus saturations
        minPc_ = MaterialLaw::endPointPc(materialParams);
        pn_ = problem.nonWettingReferencePressure();
        //porosity
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, ParentType::numComponents());
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, fluidSystemPhaseIdx);

        const int compIIdx = fluidSystemPhaseIdx;
        for (unsigned int compJIdx = 0; compJIdx < ParentType::numComponents(); ++compJIdx)
            if(compIIdx != compJIdx)
                setDiffusionCoefficient_(compJIdx,
                                         FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                 paramCache,
                                                                                 fluidSystemPhaseIdx,
                                                                                 compIIdx,
                                                                                 compJIdx));
    }

    /*!
     * \brief Fill the fluid state according to the primary variables.
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

        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const auto& priVars = elemSol[scv.localDofIndex()];

        // set the wetting pressure
        fluidState.setPressure(fluidSystemPhaseIdx, priVars[Indices::pressureIdx]);

        // compute the capillary pressure to compute the saturation
        // make sure that we the capillary pressure is not smaller than the minimum pc
        // this would possibly return unphysical values from regularized material laws
        using std::max;
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const Scalar pc = max(MaterialLaw::endPointPc(materialParams),
                              problem.nonWettingReferencePressure() - fluidState.pressure(fluidSystemPhaseIdx));
        const Scalar sw = MaterialLaw::sw(materialParams, pc);
        fluidState.setSaturation(fluidSystemPhaseIdx, sw);

        // set the mole/mass fractions
        if(useMoles)
        {
            Scalar sumSecondaryFractions = 0.0;
            for (int compIdx = 1; compIdx < ParentType::numComponents(); ++compIdx)
            {
                fluidState.setMoleFraction(fluidSystemPhaseIdx, compIdx, priVars[compIdx]);
                sumSecondaryFractions += priVars[compIdx];
            }
            fluidState.setMoleFraction(fluidSystemPhaseIdx, 0, 1.0 - sumSecondaryFractions);
        }
        else
        {
            for (int compIdx = 1; compIdx < ParentType::numComponents(); ++compIdx)
                fluidState.setMassFraction(fluidSystemPhaseIdx, compIdx, priVars[compIdx]);
        }

        // density and viscosity
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        fluidState.setDensity(fluidSystemPhaseIdx, FluidSystem::density(fluidState, paramCache, fluidSystemPhaseIdx));
        fluidState.setViscosity(fluidSystemPhaseIdx, FluidSystem::viscosity(fluidState, paramCache, fluidSystemPhaseIdx));

        // compute and set the enthalpy
        fluidState.setEnthalpy(fluidSystemPhaseIdx, EnergyVolVars::enthalpy(fluidState, paramCache, fluidSystemPhaseIdx));
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Return the temperature
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
    Scalar saturation(const int phaseIdx = fluidSystemPhaseIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? fluidState_.saturation(fluidSystemPhaseIdx) : 1.0-fluidState_.saturation(fluidSystemPhaseIdx); }

    /*!
     * \brief Returns the average mass density \f$\mathrm{[kg/m^3]}\f$ of a given
     *        fluid phase within the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar density(const int phaseIdx = fluidSystemPhaseIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? fluidState_.density(phaseIdx) : 0.0; }

    /*!
     * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * For the non-wetting phase (i.e. the gas phase), we assume
     * infinite mobility, which implies that the non-wetting phase
     * pressure is equal to the finite volume's reference pressure
     * defined by the problem.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar pressure(const int phaseIdx = fluidSystemPhaseIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? fluidState_.pressure(phaseIdx) : pn_; }

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
    Scalar mobility(const int phaseIdx = fluidSystemPhaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     * \note The non-wetting phase is infinitely mobile
     */
    Scalar viscosity(const int phaseIdx = fluidSystemPhaseIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? fluidState_.viscosity(fluidSystemPhaseIdx) : 0.0; }

    /*!
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar relativePermeability(const int phaseIdx = fluidSystemPhaseIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? relativePermeabilityWetting_ : 1.0; }

    /*!
     * \brief Returns the effective capillary pressure \f$\mathrm{[Pa]}\f$ within the
     *        control volume.
     *
     * The capillary pressure is defined as the difference in
     * pressures of the non-wetting and the wetting phase, i.e.
     * \f[ p_c = p_n - p_w \f]
     *
     * \note Capillary pressures are always larger than the entry pressure
     *       This regularization doesn't affect the residual in which pc is not needed.
     */
    Scalar capillaryPressure() const
    {
        using std::max;
        return max(minPc_, pn_ - fluidState_.pressure(fluidSystemPhaseIdx));
    }

    /*!
     * \brief Returns the pressureHead \f$\mathrm{[cm]}\f$ of a given phase within
     *        the control volume.
     *
     * For the non-wetting phase (i.e. the gas phase), we assume
     * infinite mobility, which implies that the non-wetting phase
     * pressure is equal to the finite volume's reference pressure
     * defined by the problem.
     *
     * \param phaseIdx The index of the fluid phase
     * \note this function is here as a convenience to the user to not have to
     *       manually do a conversion. It is not correct if the density is not constant
     *       or the gravity different
     */
    Scalar pressureHead(const int phaseIdx = fluidSystemPhaseIdx) const
    { return 100.0 *(pressure(phaseIdx) - pn_)/density(phaseIdx)/9.81; }

    /*!
     * \brief Returns the water content
     *        fluid phase within the finite volume.
     *
     * The water content is defined as the fraction of
     * the saturation devided by the porosity

     * \param phaseIdx The index of the fluid phase
     * \note this function is here as a convenience to the user to not have to
     *       manually do a conversion.
     */
    Scalar waterContent(const int phaseIdx = fluidSystemPhaseIdx) const
    { return saturation(phaseIdx) * solidState_.porosity(); }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarDensity(const int phaseIdx = fluidSystemPhaseIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? this->fluidState_.molarDensity(phaseIdx) : 0.0; }

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? this->fluidState_.moleFraction(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar massFraction(const int phaseIdx, const int compIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? this->fluidState_.massFraction(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarity(const int phaseIdx, const int compIdx) const
    { return phaseIdx == fluidSystemPhaseIdx ? this->fluidState_.molarity(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component
     */
    Scalar diffusionCoefficient(const int phaseIdx, const int compIdx) const
    { return diffCoefficient_[compIdx-1]; }

protected:
    FluidState fluidState_; //!< the fluid state

private:
    /*!
     * \brief TODO docme!
     *
     * \param d TODO docme!
     * \param compIdx The index of the component
     */
    void setDiffusionCoefficient_(int compIdx, Scalar d)
    { diffCoefficient_[compIdx-1] = d; }

    std::array<Scalar, ParentType::numComponents()-1> diffCoefficient_;

    Scalar relativePermeabilityWetting_; //!< the relative permeability of the wetting phase
    SolidState solidState_;
    PermeabilityType permeability_; //!< the instrinsic permeability
    Scalar pn_; //!< the reference non-wetting pressure
    Scalar minPc_; //!< the minimum capillary pressure (entry pressure)
};

} // end namespace Dumux

#endif
