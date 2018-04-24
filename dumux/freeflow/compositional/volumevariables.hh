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
 * \ingroup FreeflowNCModel
 *
 * \copydoc Dumux::FreeflowNCVolumeVariables
 */
#ifndef DUMUX_FREEFLOW_NC_VOLUMEVARIABLES_HH
#define DUMUX_FREEFLOW_NC_VOLUMEVARIABLES_HH

#include <dune/common/exceptions.hh>

#include <dumux/freeflow/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowNCModel
 * \brief Volume variables for the single-phase, multi-component free-flow model.
 */
template <class Traits>
class FreeflowNCVolumeVariables : public FreeFlowVolumeVariables< Traits, FreeflowNCVolumeVariables<Traits> >
{
    using ThisType = FreeflowNCVolumeVariables<Traits>;
    using ParentType = FreeFlowVolumeVariables<Traits, ThisType>;

    using Scalar = typename Traits::PrimaryVariables::value_type;
    using Indices = typename Traits::ModelTraits::Indices;

    static constexpr int fluidSystemPhaseIdx = Indices::fluidSystemPhaseIdx;
    static constexpr int numComponents = Traits::ModelTraits::numComponents();

    static constexpr bool useMoles = Traits::ModelTraits::useMoles();

public:
    //! export the underlying fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! export the fluid state type
    using FluidState = typename Traits::FluidState;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);
        for (unsigned int compIIdx = 0; compIIdx < numComponents; ++compIIdx)
        {
            for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            {
                // binary diffusion coefficients
                if(compIIdx != compJIdx)
                {
                    diffCoefficient_[compIIdx][compJIdx]
                        = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                  paramCache,
                                                                  fluidSystemPhaseIdx,
                                                                  compIIdx,
                                                                  compJIdx);
                }
            }
        }
    }

    /*!
     * \brief Update the fluid state
     */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    static void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        fluidState.setTemperature(ParentType::temperature(elemSol, problem, element, scv));
        fluidState.setPressure(fluidSystemPhaseIdx, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(fluidSystemPhaseIdx, 1.0);

        Scalar sumFracMinorComp = 0.0;

        for(int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if(compIdx == Indices::mainCompIdx)
                continue;

            int offset = Indices::mainCompIdx != 0 ? 1 : 0;

            // temporary add 1.0 to remove spurious differences in mole fractions
            // which are below the numerical accuracy
            Scalar moleOrMassFraction = elemSol[0][Indices::conti0EqIdx+compIdx+offset] + 1.0;
            moleOrMassFraction = moleOrMassFraction - 1.0;
            if(useMoles)
                fluidState.setMoleFraction(fluidSystemPhaseIdx, compIdx, moleOrMassFraction);
            else
                fluidState.setMassFraction(fluidSystemPhaseIdx, compIdx, moleOrMassFraction);
            sumFracMinorComp += moleOrMassFraction;
        }
        if(useMoles)
            fluidState.setMoleFraction(fluidSystemPhaseIdx, Indices::mainCompIdx, 1.0 - sumFracMinorComp);
        else
            fluidState.setMassFraction(fluidSystemPhaseIdx, Indices::mainCompIdx, 1.0 - sumFracMinorComp);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, fluidSystemPhaseIdx);

        Scalar value = FluidSystem::density(fluidState, paramCache, fluidSystemPhaseIdx);
        fluidState.setDensity(fluidSystemPhaseIdx, value);

        value = FluidSystem::viscosity(fluidState, paramCache, fluidSystemPhaseIdx);
        fluidState.setViscosity(fluidSystemPhaseIdx, value);

        // compute and set the enthalpy
        const Scalar h = ParentType::enthalpy(fluidState, paramCache);
        fluidState.setEnthalpy(fluidSystemPhaseIdx, h);
    }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     */
    Scalar pressure() const
    { return fluidState_.pressure(fluidSystemPhaseIdx); }

    /*!
     * \brief Return the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase within the
     *        control volume.
     */
    Scalar density() const
    { return fluidState_.density(fluidSystemPhaseIdx); }

    /*!
     * \brief Return temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperatures of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Return the effective dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar effectiveViscosity() const
    { return viscosity(); }

    /*!
     * \brief Return the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid within the
     *        control volume.
     */
    Scalar viscosity() const
    { return fluidState_.viscosity(fluidSystemPhaseIdx); }

    /*!
     * \brief Returns the mass fraction of a component in the phase \f$\mathrm{[-]}\f$
     *
     * \param compIdx the index of the component
     */
    Scalar massFraction(int compIdx) const
    {
        return fluidState_.massFraction(fluidSystemPhaseIdx, compIdx);
    }

    /*!
     * \brief Returns the mole fraction of a component in the phase \f$\mathrm{[-]}\f$
     *
     * \param compIdx the index of the component
     */
    Scalar moleFraction(int compIdx) const
    {
        return fluidState_.moleFraction(fluidSystemPhaseIdx, compIdx);
    }

    /*!
     * \brief Returns the mass density of a given phase \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar molarDensity() const
    {
        return fluidState_.molarDensity(fluidSystemPhaseIdx);
    }

     /*!
     * \brief Returns the diffusion coefficient \f$\mathrm{[m^2/s]}\f$
     *
     * \param compIIdx the index of the component which diffusive
     * \param compJIdx the index of the component with respect to which compIIdx diffuses
     */
    Scalar diffusionCoefficient(int compIIdx, int compJIdx = fluidSystemPhaseIdx) const
    {
        if (compIIdx == compJIdx)
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient called for fluidSystemPhaseIdx = compIdx");
        return diffCoefficient_[compIIdx][compJIdx];
    }

     /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$
     *
     * \param compIIdx the index of the component which diffusive
     * \param compJIdx the index of the component with respect to which compIIdx diffuses
     */
    Scalar effectiveDiffusivity(int compIIdx, int compJIdx = fluidSystemPhaseIdx) const
    {
        return diffusionCoefficient(compIIdx, compJIdx);
    }

    /*!
     * \brief Return the fluid state of the control volume.
     */
    const FluidState& fluidState() const
    { return fluidState_; }

protected:
    FluidState fluidState_;
    std::array<std::array<Scalar, numComponents>, numComponents> diffCoefficient_;
};

} // end namespace Dumux

#endif
