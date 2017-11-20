// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 *
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase, n-component mineralization model.
 */
#ifndef DUMUX_2PNCMIN_VOLUME_VARIABLES_HH
#define DUMUX_2PNCMIN_VOLUME_VARIABLES_HH

#include <dumux/common/math.hh>

#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/discretization/volumevariables.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/constraintsolvers/miscible2pnccomposition.hh>
#include <dumux/porousmediumflow/2pnc/implicit/volumevariables.hh>

#include "properties.hh"
#include "indices.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPNCMinModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, n-component model.
 */
template <class TypeTag>
class TwoPNCMinVolumeVariables : public TwoPNCVolumeVariables<TypeTag>
{
    // base type is used for energy related quantities
    using BaseType = ImplicitVolumeVariables<TypeTag>;

    using ParentType = TwoPNCVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum
    {
        dim = GridView::dimension,
        dimWorld=GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases =  GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents),

        // formulations
        formulation = GET_PROP_VALUE(TypeTag, Formulation),
        pwsn = TwoPNCFormulation::pwsn,
        pnsw = TwoPNCFormulation::pnsw,

        // phase indices
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        // component indices
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        // phase presence enums
        nPhaseOnly = Indices::nPhaseOnly,
        wPhaseOnly = Indices::wPhaseOnly,
        bothPhases = Indices::bothPhases,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CoordScalar = typename Grid::ctype;
    using Miscible2pNCComposition = Dumux::Miscible2pNCComposition<Scalar, FluidSystem>;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FluidSystem>;

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \copydoc ImplicitVolumeVariables::update
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        // Update parent type (also completes the fluid state)
        ParentType::update(elemSol, problem, element, scv);

        Implementation::completeFluidState(elemSol, problem, element, scv, fluidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        auto&& priVars = elemSol[scv.indexInElement()];

        sumPrecipitates_ = 0.0;
        for(int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
        {
           precipitateVolumeFraction_[sPhaseIdx] = priVars[numComponents + sPhaseIdx];
           sumPrecipitates_+= precipitateVolumeFraction_[sPhaseIdx];
        }

        salinity_= 0.0;
        moleFractionSalinity_ = 0.0;
        for (int compIdx = numMajorComponents; compIdx< numComponents; compIdx++)    //sum of the mass fraction of the components
        {
            if(this->fluidState_.moleFraction(wPhaseIdx, compIdx) > 0)
            {
                salinity_+= this->fluidState_.massFraction(wPhaseIdx, compIdx);
                //TODO: we should be consistent. Assumin all salinity to be NaCl isn't consistent with
                // the calculation of moleFractionSalinity_ here!
                moleFractionSalinity_ += this->fluidState_.moleFraction(wPhaseIdx, compIdx);
            }
        }
    }

    /*!
     * \copydoc ImplicitModel::completeFluidState
     * \param isOldSol Specifies whether this is the previous solution or the current one
     */
    static void completeFluidState(const ElementSolutionVector& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)

    {
        Scalar t = BaseType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);
        const auto phasePresence = priVars.state();

        /////////////
        // set the saturations
        /////////////

        Scalar Sn;
        if (phasePresence == nPhaseOnly)
        {
            Sn = 1.0;
        }
        else if (phasePresence == wPhaseOnly)
        {
            Sn = 0.0;
        }
        else if (phasePresence == bothPhases)
        {
            if (formulation == pwsn)
                Sn = priVars[switchIdx];
            else if (formulation == pnsw)
                Sn = 1.0 - priVars[switchIdx];
            else
                DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");
        }

        fluidState.setSaturation(nPhaseIdx, Sn);
        fluidState.setSaturation(wPhaseIdx, 1.0 - Sn);

        /////////////
        // set the pressures of the fluid phases
        /////////////

        // calculate capillary pressure
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        auto pc = MaterialLaw::pc(materialParams, 1 - Sn);

        // extract the pressures
        if (formulation == pwsn)
        {
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx]);
            if (priVars[pressureIdx] + pc < 0.0)
                 DUNE_THROW(Dumux::NumericalProblem, "Capillary pressure is too low");
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx] + pc);
        }
        else if (formulation == pnsw)
        {
            fluidState.setPressure(nPhaseIdx, priVars[pressureIdx]);
            // Here we check for (p_g - pc) in order to ensure that (p_l > 0)
            if (priVars[pressureIdx] - pc < 0.0)
            {
                std::cout<< " The gas pressure is p_g = "<< priVars[pressureIdx]<<", the capillary pressure p_c = "<< pc << std::endl;
                DUNE_THROW(Dumux::NumericalProblem, "Capillary pressure is too high");
            }
            fluidState.setPressure(wPhaseIdx, priVars[pressureIdx] - pc);
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "Formulation: " << formulation << " is invalid.");
        }

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::ParameterCache paramCache;

        // now comes the tricky part: calculate phase composition
        if (phasePresence == bothPhases)
        {
            // both phases are present, phase composition results from
            // the nonwetting <-> wetting equilibrium. This is
            // the job of the "MiscibleMultiPhaseComposition"
            // constraint solver

            // set the known mole fractions in the fluidState so that they
            // can be used by the Miscible2pNcComposition constraint solver
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, priVars[compIdx]);
            }

            Miscible2pNCComposition::solve(fluidState,
                                           paramCache,
                                           wPhaseIdx,  //known phaseIdx
                                           /*setViscosity=*/true,
                                           /*setInternalEnergy=*/false);
        }
        else if (phasePresence == nPhaseOnly)
        {
            Dune::FieldVector<Scalar, numComponents> moleFrac;
            Dune::FieldVector<Scalar, numComponents> fugCoeffL;
            Dune::FieldVector<Scalar, numComponents> fugCoeffG;

            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fugCoeffL[compIdx] = FluidSystem::fugacityCoefficient(fluidState,
                                                                      paramCache,
                                                                      wPhaseIdx,
                                                                      compIdx);
                fugCoeffG[compIdx] = FluidSystem::fugacityCoefficient(fluidState,
                                                                      paramCache,
                                                                      nPhaseIdx,
                                                                      compIdx);
            }
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                moleFrac[compIdx] = (priVars[compIdx]*fugCoeffL[compIdx]*fluidState.pressure(wPhaseIdx))
                                    /(fugCoeffG[compIdx]*fluidState.pressure(nPhaseIdx));

            moleFrac[wCompIdx] =  priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                sumMoleFracOtherComponents+=moleFrac[compIdx];

            sumMoleFracOtherComponents += moleFrac[wCompIdx];
            moleFrac[nCompIdx] = 1 - sumMoleFracOtherComponents;

            // Set fluid state mole fractions
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(nPhaseIdx, compIdx, moleFrac[compIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                                    paramCache,
                                                    nPhaseIdx,
//                                                     nPhaseOnly,
                                                    /*setViscosity=*/true,
                                                    /*setEnthalpy=*/false);

        }
        else if (phasePresence == wPhaseOnly)
        {
            // only the wetting phase is present, i.e. wetting
            // composition is stored explicitly.
            // extract _mass_ fractions in the nonwetting phase
            Dune::FieldVector<Scalar, numComponents> moleFrac;

            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
            }
            moleFrac[nCompIdx] = priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = 0;
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                    sumMoleFracOtherComponents+=moleFrac[compIdx];
            }
            sumMoleFracOtherComponents += moleFrac[nCompIdx];
            moleFrac[wCompIdx] = 1 -sumMoleFracOtherComponents;

            // convert mass to mole fractions and set the fluid state
            for (int compIdx=0; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, moleFrac[compIdx]);
            }

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                                    paramCache,
                                                    wPhaseIdx,
//                                                     wPhaseOnly,
                                                    /*setViscosity=*/true,
                                                    /*setEnthalpy=*/false);
        }
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            Scalar h = BaseType::enthalpy(fluidState, paramCache, phaseIdx);

            fluidState.setDensity(phaseIdx, rho);
            fluidState.setViscosity(phaseIdx, mu);
            fluidState.setEnthalpy(phaseIdx, h);
        }
    }
    /*!
     * \brief Returns the volume fraction of the precipitate (solid phase)
     * for the given phaseIdx
     *
     * \param phaseIdx the index of the solid phase
     */
    Scalar precipitateVolumeFraction(int phaseIdx) const
    { return precipitateVolumeFraction_[phaseIdx - numPhases]; }

    /*!
     * \brief Returns the density of the phase for all fluid and solid phases
     *
     * \param phaseIdx the index of the fluid phase
     */
    Scalar density(int phaseIdx) const
    {
        if (phaseIdx < numPhases)
            return this->fluidState_.density(phaseIdx);
        else
            return FluidSystem::precipitateDensity(phaseIdx);
    }
    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx < numPhases)
            return this->fluidState_.molarDensity(phaseIdx);
        else
            return FluidSystem::precipitateMolarDensity(phaseIdx);
    }

    /*!
     * \brief Returns the molality of a component in the phase
     *
     * \param phaseIdx the index of the fluid phase
     * \param compIdx the index of the component
     * \f$\mathrm{molality}=\frac{n_\mathrm{component}}{m_\mathrm{solvent}}
     * =\frac{n_\mathrm{component}}{n_\mathrm{solvent}*M_\mathrm{solvent}}\f$
     * compIdx of the main component (solvent) in the
     * phase is equal to the phaseIdx
     */
    Scalar molality(int phaseIdx, int compIdx) const // [moles/Kg]
    {
        return this->fluidState_.moleFraction(phaseIdx, compIdx)
                  /(this->fluidState_.moleFraction(phaseIdx, phaseIdx)
                  * FluidSystem::molarMass(phaseIdx));
    }

protected:

    Scalar precipitateVolumeFraction_[numSPhases];
    Scalar sumPrecipitates_;
    FluidState fluidState_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace

#endif
