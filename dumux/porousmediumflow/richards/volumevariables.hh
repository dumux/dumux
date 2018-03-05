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
 * \ingroup RichardsModel
 * \brief Volume averaged quantities required by the Richards model.
 */
#ifndef DUMUX_RICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_VOLUME_VARIABLES_HH

#include <cassert>

#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/material/idealgas.hh>
#include <dumux/material/constants.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Volume averaged quantities required by the Richards model.
 *
 * This contains the quantities which are are constant within a finite
 * volume in the Richards model.
 */
template <class TypeTag>
class RichardsVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PermeabilityType = typename SpatialParams::PermeabilityType;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum{
         pressureIdx = Indices::pressureIdx,
         switchIdx = Indices::switchIdx,
         wPhaseIdx = Indices::wPhaseIdx,
         nPhaseIdx = Indices::nPhaseIdx,
         wCompIdx = FluidSystem::wCompIdx,
         nCompIdx = FluidSystem::nCompIdx
    };

    // present phases
    enum
    {
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool enableWaterDiffusionInAir
        = GET_PROP_VALUE(TypeTag, EnableWaterDiffusionInAir);
    static constexpr bool useKelvinVaporPressure
        = GET_PROP_VALUE(TypeTag, UseKelvinEquation);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        static_assert(!(!enableWaterDiffusionInAir && useKelvinVaporPressure),
          "Kevin vapor presssure only makes sense if water in air is considered!");

        ParentType::update(elemSol, problem, element, scv);
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);
        const auto phasePresence = priVars.state();

        // precompute the minimum capillary pressure (entry pressure)
        // needed to make sure we don't compute unphysical capillary pressures and thus saturations
        minPc_ = MaterialLaw::endPointPc(materialParams);

        if (phasePresence == nPhaseOnly)
        {
            moleFraction_[wPhaseIdx] = 1.0;
            moleFraction_[nPhaseIdx] = priVars[switchIdx];

            fluidState_.setSaturation(wPhaseIdx, 0.0);
            fluidState_.setSaturation(nPhaseIdx, 1.0);

            Scalar t = ParentType::temperature(elemSol, problem, element, scv);
            fluidState_.setTemperature(t);

            // get pc for sw = 0.0
            const Scalar pc = MaterialLaw::pc(materialParams, 0.0);

            // set the wetting pressure
            fluidState_.setPressure(wPhaseIdx, problem.nonWettingReferencePressure() - pc);
            fluidState_.setPressure(nPhaseIdx, problem.nonWettingReferencePressure());

            // set molar densities
            molarDensity_[wPhaseIdx] = FluidSystem::H2O::liquidDensity(temperature(), pressure(wPhaseIdx))/FluidSystem::H2O::molarMass();
            molarDensity_[nPhaseIdx] = IdealGas<Scalar>::molarDensity(temperature(), problem.nonWettingReferencePressure());

            // density and viscosity
            typename FluidSystem::ParameterCache paramCache;
            paramCache.updateAll(fluidState_);
            fluidState_.setDensity(wPhaseIdx, FluidSystem::density(fluidState_, paramCache, wPhaseIdx));
            fluidState_.setDensity(nPhaseIdx, FluidSystem::density(fluidState_, paramCache, nPhaseIdx));
            fluidState_.setViscosity(wPhaseIdx, FluidSystem::viscosity(fluidState_, paramCache, wPhaseIdx));

            // compute and set the enthalpy
            fluidState_.setEnthalpy(wPhaseIdx, Implementation::enthalpy(fluidState_, paramCache, wPhaseIdx));
            fluidState_.setEnthalpy(nPhaseIdx, Implementation::enthalpy(fluidState_, paramCache, nPhaseIdx));
        }
        else if (phasePresence == bothPhases)
        {
            Implementation::completeFluidState(elemSol, problem, element, scv, fluidState_);

            // if we want to account for diffusion in the air phase
            // use Raoult to compute the water mole fraction in air
            if (enableWaterDiffusionInAir)
            {
                molarDensity_[wPhaseIdx] = FluidSystem::H2O::liquidDensity(temperature(), pressure(wPhaseIdx))/FluidSystem::H2O::molarMass();
                molarDensity_[nPhaseIdx] = IdealGas<Scalar>::molarDensity(temperature(), problem.nonWettingReferencePressure());
                moleFraction_[wPhaseIdx] = 1.0;

                moleFraction_[nPhaseIdx] = FluidSystem::H2O::vaporPressure(temperature()) / problem.nonWettingReferencePressure();
                if (useKelvinVaporPressure)
                {
                    using std::exp;
                    moleFraction_[nPhaseIdx] *= exp(-capillaryPressure() * FluidSystem::H2O::molarMass()/density(wPhaseIdx)
                                                     / Constants<Scalar>::R / temperature());
                }

                // binary diffusion coefficients
                typename FluidSystem::ParameterCache paramCache;
                paramCache.updateAll(fluidState_);
                diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(fluidState_, paramCache, nPhaseIdx, wCompIdx, nCompIdx);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            Implementation::completeFluidState(elemSol, problem, element, scv, fluidState_);
        }

        //////////
        // specify the other parameters
        //////////
        relativePermeabilityWetting_ = MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx));
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
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
    template<class ElementSolution>
    static void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);

        // set the wetting pressure
        using std::max;
        Scalar minPc = MaterialLaw::pc(materialParams, 1.0);
        fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
        fluidState.setPressure(nPhaseIdx, max(problem.nonWettingReferencePressure(), fluidState.pressure(wPhaseIdx) + minPc));

        // compute the capillary pressure to compute the saturation
        // make sure that we the capillary pressure is not smaller than the minimum pc
        // this would possibly return unphysical values from regularized material laws
        using std::max;
        const Scalar pc = max(MaterialLaw::endPointPc(materialParams),
                              problem.nonWettingReferencePressure() - fluidState.pressure(wPhaseIdx));
        const Scalar sw = MaterialLaw::sw(materialParams, pc);
        fluidState.setSaturation(wPhaseIdx, sw);
        fluidState.setSaturation(nPhaseIdx, 1.0-sw);

        // density and viscosity
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        fluidState.setDensity(wPhaseIdx, FluidSystem::density(fluidState, paramCache, wPhaseIdx));
        fluidState.setDensity(nPhaseIdx, FluidSystem::density(fluidState, paramCache, nPhaseIdx));
        fluidState.setViscosity(wPhaseIdx, FluidSystem::viscosity(fluidState, paramCache, wPhaseIdx));

        // compute and set the enthalpy
        fluidState.setEnthalpy(wPhaseIdx, Implementation::enthalpy(fluidState, paramCache, wPhaseIdx));
        fluidState.setEnthalpy(nPhaseIdx, Implementation::enthalpy(fluidState, paramCache, nPhaseIdx));
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

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
    { return porosity_; }

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
    Scalar saturation(const int phaseIdx = wPhaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the average mass density \f$\mathrm{[kg/m^3]}\f$ of a given
     *        fluid phase within the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar density(const int phaseIdx = wPhaseIdx) const
    { return fluidState_.density(phaseIdx); }

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
    Scalar pressure(const int phaseIdx = wPhaseIdx) const
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
    Scalar mobility(const int phaseIdx = wPhaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Returns the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     * \note The non-wetting phase is infinitely mobile
     */
    Scalar viscosity(const int phaseIdx = wPhaseIdx) const
    { return phaseIdx == wPhaseIdx ? fluidState_.viscosity(wPhaseIdx) : 0.0; }

    /*!
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar relativePermeability(const int phaseIdx = wPhaseIdx) const
    { return phaseIdx == wPhaseIdx ? relativePermeabilityWetting_ : 1.0; }

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
        return max(minPc_, pressure(nPhaseIdx) - pressure(wPhaseIdx));
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
    Scalar pressureHead(const int phaseIdx = wPhaseIdx) const
    { return 100.0 *(pressure(phaseIdx) - pressure(nPhaseIdx))/density(phaseIdx)/9.81; }

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
    Scalar waterContent(const int phaseIdx = wPhaseIdx) const
    { return saturation(phaseIdx) * porosity_; }

    /*!
     * \brief Returns the mole fraction of a given component in a
     *        given phase within the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    {
        assert(enableWaterDiffusionInAir);
        if (compIdx != wCompIdx)
            DUNE_THROW(Dune::InvalidStateException, "There is only one component for Richards!");
        return moleFraction_[phaseIdx];
    }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume in \f$[mol/m^3]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(const int phaseIdx) const
    {
        assert(enableWaterDiffusionInAir);
        return molarDensity_[phaseIdx];
    }

    /*!
     * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
     *
     * \param phaseIdx The phase index
     * \param compIdx The component index
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        assert(enableWaterDiffusionInAir && phaseIdx == nPhaseIdx && compIdx == wCompIdx);
        return diffCoeff_;
    }

protected:
    FluidState fluidState_; //!< the fluid state
    Scalar relativePermeabilityWetting_; //!< the relative permeability of the wetting phase
    Scalar porosity_; //!< the porosity
    PermeabilityType permeability_; //!< the instrinsic permeability
    Scalar minPc_; //!< the minimum capillary pressure (entry pressure)
    Scalar moleFraction_[numPhases]; //!< The water mole fractions in water and air
    Scalar molarDensity_[numPhases]; //!< The molar density of water and air
    Scalar diffCoeff_; //!< The binary diffusion coefficient of water in air
};

} // end namespace Dumux

#endif
