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

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/volumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsNCModel
 * \brief  Contains the quantities which are constant within a
 *        finite volume in the Richards, n-component model.
 */
template <class TypeTag>
class RichardsBaseVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
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
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum{
         pressureIdx = Indices::pressureIdx,
         wPhaseIdx = Indices::wPhaseIdx
    };

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

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
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        Implementation::completeFluidState(elemSol, problem, element, scv, fluidState_);
        //////////
        // specify the other parameters
        //////////
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        relativePermeabilityWetting_ = MaterialLaw::krw(materialParams, fluidState_.saturation(wPhaseIdx));

        // precompute the minimum capillary pressure (entry pressure)
        // needed to make sure we don't compute unphysical capillary pressures and thus saturations
        minPc_ = MaterialLaw::endPointPc(materialParams);
        pn_ = problem.nonWettingReferencePressure();
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
    static void completeFluidState(const ElementSolutionVector& elemSol,
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
        fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);

        // compute the capillary pressure to compute the saturation
        // make sure that we the capillary pressure is not smaller than the minimum pc
        // this would possibly return unphysical values from regularized material laws
        using std::max;
        const Scalar pc = max(MaterialLaw::endPointPc(materialParams),
                              problem.nonWettingReferencePressure() - fluidState.pressure(wPhaseIdx));
        const Scalar sw = MaterialLaw::sw(materialParams, pc);
        fluidState.setSaturation(wPhaseIdx, sw);

        // density and viscosity
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        fluidState.setDensity(wPhaseIdx, FluidSystem::density(fluidState, paramCache, wPhaseIdx));
        fluidState.setViscosity(wPhaseIdx, FluidSystem::viscosity(fluidState, paramCache, wPhaseIdx));

        // compute and set the enthalpy
        fluidState.setEnthalpy(wPhaseIdx, Implementation::enthalpy(fluidState, paramCache, wPhaseIdx));
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
    { return phaseIdx == wPhaseIdx ? fluidState_.saturation(wPhaseIdx) : 1.0-fluidState_.saturation(wPhaseIdx); }

    /*!
     * \brief Returns the average mass density \f$\mathrm{[kg/m^3]}\f$ of a given
     *        fluid phase within the control volume.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar density(const int phaseIdx = wPhaseIdx) const
    { return phaseIdx == wPhaseIdx ? fluidState_.density(phaseIdx) : 0.0; }

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
    { return phaseIdx == wPhaseIdx ? fluidState_.pressure(phaseIdx) : pn_; }

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
        return max(minPc_, pn_ - fluidState_.pressure(wPhaseIdx));
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
    Scalar waterContent(const int phaseIdx = wPhaseIdx) const
    { return saturation(phaseIdx) * porosity_; }

protected:
    FluidState fluidState_; //!< the fluid state
    Scalar relativePermeabilityWetting_; //!< the relative permeability of the wetting phase
    Scalar porosity_; //!< the porosity
    PermeabilityType permeability_; //!< the instrinsic permeability
    Scalar pn_; //!< the reference non-wetting pressure
    Scalar minPc_; //!< the minimum capillary pressure (entry pressure)
};

/*!
 * \ingroup RichardsNCModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the Richards, n-component model.
 */
template <class TypeTag>
class RichardsNCVolumeVariables : public RichardsBaseVolumeVariables<TypeTag>
{
    using ParentType = RichardsBaseVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        pressureIdx = Indices::pressureIdx
    };

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static const int dimWorld = GridView::dimensionworld;
    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

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
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        //calculate all secondary variables from the primary variables and store results in fluidstate
        Implementation::completeFluidState(elemSol, problem, element, scv, this->fluidState_);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(this->fluidState_, wPhaseIdx);

        const int compIIdx = wPhaseIdx;
        for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            if(compIIdx != compJIdx)
                setDiffusionCoefficient_(compJIdx,
                                         FluidSystem::binaryDiffusionCoefficient(this->fluidState_,
                                                                                 paramCache,
                                                                                 wPhaseIdx,
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
     * \param elemSol A vector containing all primary variables connected to the element.
     * \param problem The problem at hand.
     * \param element The current element.
     * \param scv The subcontrol volume.
     * \param fluidState The fluid state to fill.
     */
    static void completeFluidState(const ElementSolutionVector &elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume &scv,
                                   FluidState& fluidState)
    {
        ParentType::completeFluidState(elemSol, problem, element, scv, fluidState);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);

        // set the mole/mass fractions
        if(useMoles)
        {
            Scalar sumSecondaryFractions = 0.0;
            for (int compIdx = 1; compIdx < numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, priVars[compIdx]);
                sumSecondaryFractions += priVars[compIdx];
            }
            fluidState.setMoleFraction(wPhaseIdx, 0, 1.0 - sumSecondaryFractions);
        }
        else
        {
            for (int compIdx = 1; compIdx < numComponents; ++compIdx)
                fluidState.setMassFraction(wPhaseIdx, compIdx, priVars[compIdx]);
        }
    }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarDensity(const int phaseIdx = wPhaseIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.molarDensity(phaseIdx) : 0.0; }

    /*!
     * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar moleFraction(const int phaseIdx, const int compIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.moleFraction(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return mass fraction \f$\mathrm{[kg/kg]}\f$ of a component in the phase.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar massFraction(const int phaseIdx, const int compIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.massFraction(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return concentration \f$\mathrm{[mol/m^3]}\f$  of a component in the phase.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarity(const int phaseIdx, const int compIdx) const
    { return phaseIdx == wPhaseIdx ? this->fluidState_.molarity(phaseIdx, compIdx) : 0.0; }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     *
     * \param phaseIdx The index of the phase.
     * \param compIdx The index of the component
     */
    Scalar diffusionCoefficient(const int phaseIdx, const int compIdx) const
    {
        assert(phaseIdx == wPhaseIdx);
        assert(compIdx > wPhaseIdx);
        return diffCoefficient_[compIdx-1];
    }

    /*!
     * \brief Returns the dispersivity of the fluid's streamlines.
     */
    const GlobalPosition &dispersivity() const
    { return dispersivity_; }

private:
    /*!
     * \brief TODO docme!
     *
     * \param d TODO docme!
     * \param compIdx The index of the component
     */
    void setDiffusionCoefficient_(int compIdx, Scalar d)
    {
        assert(compIdx > wPhaseIdx);
        diffCoefficient_[compIdx-1] = d;
    }

    std::array<Scalar, numComponents-1> diffCoefficient_;
    GlobalPosition dispersivity_;
};

} // end namespace Dumux

#endif
