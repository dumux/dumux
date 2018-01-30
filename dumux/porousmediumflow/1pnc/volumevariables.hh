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
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/volumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup OnePNCModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase, n-component model.
 *
 * \note The return functions for the fluid state variables always forward to the actual
 *       fluid state using the phaseIdx from the DuMuX property system. Furthermore, the
 *       default value is not used, but is only here to enable calling these functions
 *       without handing in a phase index (as in a single-phasic context there is only one phase).
 *       This way one can use two-phase fluid systems for this one-phasic flow and transport
 *       model by specifying which phase is present through the DuMuX property system.
 */
template <class TypeTag>
class OnePNCVolumeVariables : public PorousMediumFlowVolumeVariables<TypeTag>
{
    using ParentType = PorousMediumFlowVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PermeabilityType = typename SpatialParams::PermeabilityType;

    enum
    {
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        phaseIdx = Indices::phaseIdx,
        phaseCompIdx = Indices::phaseCompIdx,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        firstMoleFracIdx = Indices::firstMoleFracIdx,

    };

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

        completeFluidState(elemSol, problem, element, scv, fluidState_);

        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, phaseIdx);

        int compIIdx = phaseCompIdx;

        for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        {
            diffCoeff_[compJIdx] = 0.0;
            if(compIIdx!= compJIdx)
                {
                diffCoeff_[compJIdx] = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                             paramCache,
                                                             phaseIdx,
                                                             compIIdx,
                                                             compJIdx);
                }
        }

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
    static void completeFluidState(const ElementSolutionVector &elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume &scv,
                                   FluidState& fluidState)

    {
        Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);
        fluidState.setSaturation(phaseIdx, 1.);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);
        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        // calculate the phase composition
        Dune::FieldVector<Scalar, numComponents> moleFrac;

        Scalar sumMoleFracNotWater = 0;

        for (int compIdx=firstMoleFracIdx; compIdx<numComponents; ++compIdx){
             moleFrac[compIdx] = priVars[compIdx];
             sumMoleFracNotWater +=moleFrac[compIdx];
            }
        moleFrac[0] = 1- sumMoleFracNotWater;

        // Set fluid state mole fractions
        for (int compIdx=0; compIdx<numComponents; ++compIdx)
        {
            fluidState.setMoleFraction(phaseIdx, compIdx, moleFrac[compIdx]);
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
        Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);

        fluidState.setDensity(phaseIdx, rho);
        fluidState.setViscosity(phaseIdx, mu);

        // compute and set the enthalpy
        Scalar h = Implementation::enthalpy(fluidState, paramCache, phaseIdx);
        fluidState.setEnthalpy(phaseIdx, h);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models. We always forward to the fluid state with the
     *       phaseIdx property (see class description).
     */
    Scalar density(int pIdx = phaseIdx) const
    {
        assert(pIdx == phaseIdx);
        return fluidState_.density(phaseIdx);
    }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models. We always forward to the fluid state with the
     *       phaseIdx property (see class description).
     */
    Scalar molarDensity(int pIdx = phaseIdx) const
    {
        assert(pIdx == phaseIdx);
        return fluidState_.molarDensity(phaseIdx);
    }

    /*!
     * \brief Return the saturation
     *
     * This method is here for compatibility reasons with other models. The saturation
     * is always 1.0 in a one-phasic context.
     */
    Scalar saturation(int pIdx = phaseIdx) const
    { return 1.0; }

     /*!
      * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      *
      * \note the phase index passed to this function is for compatibility reasons
      *       with multiphasic models. We always forward to the fluid state with the
      *       phaseIdx property (see class description).
      */
     Scalar moleFraction(int pIdx, int compIdx) const
     {
         // make sure this is only called with admissible indices
         assert(pIdx == phaseIdx);
         assert(compIdx < numComponents);
         return fluidState_.moleFraction(phaseIdx, compIdx);
     }

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      *
      * \note the phase index passed to this function is for compatibility reasons
      *       with multiphasic models. We always forward to the fluid state with the
      *       phaseIdx property (see class description).
      */
     Scalar massFraction(int pIdx, int compIdx) const
     {
         // make sure this is only called with admissible indices
         assert(pIdx == phaseIdx);
         assert(compIdx < numComponents);
         return fluidState_.massFraction(phaseIdx, compIdx);
     }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models. We always forward to the fluid state with the
     *       phaseIdx property (see class description).
     */
    Scalar pressure(int pIdx = phaseIdx) const
    {
        assert(pIdx == phaseIdx);
        return fluidState_.pressure(phaseIdx);
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
     *       with multiphasic models. We always forward to the fluid state with the
     *       phaseIdx property (see class description).
     */
    Scalar mobility(int pIdx = phaseIdx) const
    {
        assert(pIdx == phaseIdx);
        return 1.0/fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     *
     * \note the phase index passed to this function is for compatibility reasons
     *       with multiphasic models. We always forward to the fluid state with the
     *       phaseIdx property (see class description).
     */
    Scalar viscosity(int pIdx = phaseIdx) const
    {
        assert(pIdx == phaseIdx);
        return fluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffusionCoefficient(int pIdx, int compIdx) const
    {
        assert(pIdx == phaseIdx);
        assert(compIdx < numComponents);
        return diffCoeff_[compIdx];
    }

    /*!
     * \brief Returns the molarity of a component in the phase
     *
     * \param compIdx the index of the component
     */
    Scalar molarity(int compIdx) const // [moles/m^3]
    {
        assert(compIdx < numComponents);
        return fluidState_.molarity(phaseIdx, compIdx);
    }

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param compIdx the index of the component
      */
     Scalar massFraction(int compIdx) const
     {
         assert(compIdx < numComponents);
         return this->fluidState_.massFraction(phaseIdx, compIdx);
     }

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

protected:

    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               int scvIdx)
    {
         return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }


    Scalar porosity_;        //!< Effective porosity within the control volume
    PermeabilityType permeability_;
    Scalar density_;
    FluidState fluidState_;
    Dune::FieldVector<Scalar, numComponents> diffCoeff_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

};

} // end namespace

#endif
