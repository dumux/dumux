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
 * \ingroup OnePNCModel
 * \brief Quantities required by the single-phase, n-component box
 *        model defined on a vertex.
 */
#ifndef DUMUX_1PNC_VOLUME_VARIABLES_HH
#define DUMUX_1PNC_VOLUME_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/porousmediumflow/volumevariables.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase, n-component model.
 *
 * \note The default value for the phase index given in the fluid property interfaces is not used,
 *       but is only here to enable calling these functions without handing in a phase index
 *       (as in a single-phasic context there is only one phase).
 */
template <class Traits>
class OnePNCVolumeVariables
: public PorousMediumFlowVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, OnePNCVolumeVariables<Traits> >
{
    using ParentType = PorousMediumFlowVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, OnePNCVolumeVariables<Traits> >;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using Idx = typename Traits::ModelTraits::Indices;
    static constexpr int numFluidComps = ParentType::numComponents();

    enum
    {
        // pressure primary variable index
        pressureIdx = Idx::pressureIdx
    };

public:
    //! export fluid state type
    using FluidState = typename Traits::FluidState;
    //! export fluid system type
    using FluidSystem = typename Traits::FluidSystem;
    //! export indices
    using Indices = typename Traits::ModelTraits::Indices;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    using SolidSystem = typename Traits::SolidSystem;

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
                const Scv &scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        // calculate the remaining quantities
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, 0);

        diffCoeff_[0] = 0.0; // the main component with itself doesn't have a binary diffusion coefficient
        for (unsigned int compJIdx = 1; compJIdx < numFluidComps; ++compJIdx)
            diffCoeff_[compJIdx] = FluidSystem::binaryDiffusionCoefficient(fluidState_, paramCache, 0, 0, compJIdx);
    }

    /*!
     * \brief Set complete fluid state
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     * \param fluidState A container with the current (physical) state of the fluid
     */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol &elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv &scv,
                            FluidState& fluidState,
                            SolidState& solidState)

    {
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);
        fluidState.setSaturation(0, 1.0);

        const auto& priVars = elemSol[scv.localDofIndex()];
        fluidState.setPressure(0, priVars[pressureIdx]);

        // Set fluid state mole fractions
        Scalar sumMoleFracNotMainComp = 0;
        for (int compIdx = 1; compIdx < numFluidComps; ++compIdx)
        {
            fluidState.setMoleFraction(0, compIdx, priVars[compIdx]);
            sumMoleFracNotMainComp += priVars[compIdx];
        }
        fluidState.setMoleFraction(0, 0, 1.0 - sumMoleFracNotMainComp);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        Scalar rho = FluidSystem::density(fluidState, paramCache, 0);
        Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, 0);
        Scalar mu = FluidSystem::viscosity(fluidState, paramCache, 0);

        fluidState.setDensity(0, rho);
        fluidState.setMolarDensity(0, rhoMolar);
        fluidState.setViscosity(0, mu);

        // compute and set the enthalpy
        Scalar h = EnergyVolVars::enthalpy(fluidState, paramCache, 0);
        fluidState.setEnthalpy(0, h);
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
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models.
     */
    Scalar density(int phaseIdx = 0) const
    {
        return fluidState_.density(0);
    }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models.
     */
    Scalar molarDensity(int phaseIdx = 0) const
    {
        return fluidState_.molarDensity(0);
    }

    /*!
     * \brief Return the saturation
     *
     * This method is here for compatibility reasons with other models. The saturation
     * is always 1.0 in a one-phasic context.
     */
    Scalar saturation(int phaseIdx = 0) const
    { return 1.0; }

     /*!
      * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      *
      * \note the phase index passed to this function is for compatibility reasons
      *       with multiphasic models.
      */
     Scalar moleFraction(int phaseIdx, int compIdx) const
     {
         // make sure this is only called with admissible indices
         assert(compIdx < numFluidComps);
         return fluidState_.moleFraction(0, compIdx);
     }

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      *
      * \note the phase index passed to this function is for compatibility reasons
      *       with multiphasic models.
      */
     Scalar massFraction(int phaseIdx, int compIdx) const
     {
         // make sure this is only called with admissible indices
         assert(compIdx < numFluidComps);
         return fluidState_.massFraction(0, compIdx);
     }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models.
     */
    Scalar pressure(int phaseIdx = 0) const
    {
        return fluidState_.pressure(0);
    }

    /*!
     * \brief Return temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the mobility \f$\mathrm{[1/(Pa s)]}\f$.
     *
     * The term mobility is usually not employed in the one phase context.
     * The method is here for compatibility reasons with other models.
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models.
     */
    Scalar mobility(int phaseIdx = 0) const
    {
        return 1.0/fluidState_.viscosity(0);
    }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models.
     */
    Scalar viscosity(int phaseIdx = 0) const
    {
        return fluidState_.viscosity(0);
    }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return solidState_.porosity(); }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        assert(compIdx < numFluidComps);
        return diffCoeff_[compIdx];
    }

    /*!
     * \brief Returns the molarity of a component in the phase
     *
     * \param compIdx the index of the component
     */
    Scalar molarity(int compIdx) const // [moles/m^3]
    {
        assert(compIdx < numFluidComps);
        return fluidState_.molarity(0, compIdx);
    }

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param compIdx the index of the component
      */
     Scalar massFraction(int compIdx) const
     {
         assert(compIdx < numFluidComps);
         return this->fluidState_.massFraction(0, compIdx);
     }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    PermeabilityType permeability_; //!< Effective permeability within the control volume
    Dune::FieldVector<Scalar, numFluidComps> diffCoeff_; //!< Binary diffusion coefficients
};

} // end namespace Dumux

#endif
