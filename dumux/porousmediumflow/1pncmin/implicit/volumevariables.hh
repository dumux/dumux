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
#ifndef DUMUX_1PNCMIN_VOLUME_VARIABLES_HH
#define DUMUX_1PNCMIN_VOLUME_VARIABLES_HH


#include <dumux/common/math.hh>
#include <dumux/implicit/model.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/porousmediumflow/1pnc/implicit/volumevariables.hh>

#include "properties.hh"
#include "indices.hh"

namespace Dumux
{

/*!
 * \ingroup OnePNCMinModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase, n-component model.
 */
template <class TypeTag>
class OnePNCMinVolumeVariables : public OnePNCVolumeVariables<TypeTag>
{
    // base type is used for energy related quantites
    using BaseType = ImplicitVolumeVariables<TypeTag>;

    using ParentType = OnePNCVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

//     typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum
    {
        dim = GridView::dimension,
        dimWorld=GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases =  GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        numSComponents = GET_PROP_VALUE(TypeTag, NumSComponents), //if there is more than 1 component in the solid phase

        phaseCompIdx = Indices::phaseCompIdx,

        // phase indices
        phaseIdx = FluidSystem::gPhaseIdx,
        cPhaseIdx = Indices::phaseIdx +1,
        hPhaseIdx = Indices::phaseIdx +2,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        firstMoleFracIdx = Indices::firstMoleFracIdx,

    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CoordScalar = typename Grid::ctype;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

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
        ParentType::update(elemSol, problem, element, scv);

//         ParentType::update(priVars, problem, element, fvGeometry, scvIdx, isOldSol);
//         completeFluidState(priVars, problem, element, fvGeometry, scvIdx, this->fluidState_, isOldSol);

        /////////////
        // calculate the remaining quantities
        /////////////

        auto&& priVars = isBox ? elemSol[scv.index()] : elemSol[0];

        // porosity evaluation
        initialPorosity_ = problem.spatialParams().porosity(element, scv);
        minimumPorosity_ = problem.spatialParams().porosityMin(element, scv);
        maximumPorosity_ = problem.spatialParams().porosityMax(element, scv);

        sumPrecipitates_ = 0.0;
        for(int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
        {
           precipitateVolumeFraction_[sPhaseIdx] = priVars[numComponents + sPhaseIdx];
           sumPrecipitates_+= precipitateVolumeFraction_[sPhaseIdx];

        }


//         energy related quantities not contained in the fluid state
        asImp_().updateEnergy_(elemSol, problem, element, scv);
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
//         Scalar t = ParentType::temperature(elemSol, problem, element, scv); //old
        Scalar t = BaseType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);

        /////////////
        // set the saturations
        /////////////

        fluidState.setSaturation(phaseIdx, 1.0 );

        /////////////
        // set the pressures of the fluid phase
        /////////////
        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv);
        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        /////////////
        // calculate the phase compositions
        /////////////

        typename FluidSystem::ParameterCache paramCache;

        Dune::FieldVector<Scalar, numComponents> moleFrac;

        Scalar sumMoleFracNotWater = 0;
        for (int compIdx=firstMoleFracIdx; compIdx<numComponents; ++compIdx)
        {
            moleFrac[compIdx] = priVars[compIdx];
            sumMoleFracNotWater+=moleFrac[compIdx];
        }

        moleFrac[0] = 1 -sumMoleFracNotWater;

//         //mole fractions for the solid phase
//         moleFrac[firstMoleFracIdx+1]= priVars[firstMoleFracIdx+1];
//         moleFrac[firstMoleFracIdx +2] = 1- moleFrac[firstMoleFracIdx+1];

        // convert mass to mole fractions and set the fluid state
        for (int compIdx=0; compIdx<numComponents; ++compIdx)
        {
            fluidState.setMoleFraction(phaseIdx, compIdx, moleFrac[compIdx]);
        }

//         // set mole fractions for the solid phase
//         fluidState.setMoleFraction(cPhaseIdx, firstMoleFracIdx+1, moleFrac[firstMoleFracIdx+1]);
//         fluidState.setMoleFraction(hPhaseIdx, firstMoleFracIdx+2, moleFrac[firstMoleFracIdx+2]);

        paramCache.updateAll(fluidState);

//         Scalar h = Implementation::enthalpy_(fluidState, paramCache, phaseIdx);
        Scalar h = BaseType::enthalpy(fluidState, paramCache, phaseIdx);
        fluidState.setEnthalpy(phaseIdx, h);
    }
    /*!
     * \brief Returns the volume fraction of the precipitate (solid phase)
     * for the given phaseIdx
     *
     * \param phaseIdx the index of the solid phase
     */
    Scalar precipitateVolumeFraction(int phaseIdx) const
    {
        return precipitateVolumeFraction_[phaseIdx - numPhases];
    }

    /*!
     * \brief Returns the inital porosity of the
     * pure, precipitate-free porous medium
     */
    Scalar initialPorosity() const
    { return initialPorosity_;}

    /*!
     * \brief Returns the inital permeability of the
     * pure, precipitate-free porous medium
     */
    Scalar initialPermeability() const
    { return initialPermeability_;}

    /*!
     * \brief Returns the factor for the reduction of the initial permeability
     * due precipitates in the porous medium
     */
    Scalar permeabilityFactor() const
    { return permeabilityFactor_; }

    /*!
     * \brief Returns the density of the phase for all fluid and solid phases
     *
     * \param phaseIdx the index of the fluid phase
     */
    Scalar density(int phaseIdx) const
    {
        if (phaseIdx <= numPhases)
            return this->fluidState_.density(phaseIdx);
        else if (phaseIdx > numPhases)
            return FluidSystem::precipitateDensity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx <  1)
            return this->fluidState_.molarDensity(phaseIdx);
        else if (phaseIdx >= 1){
            /*Attention: sPhaseIdx of the fluidsystem and the model can be different.*/
//             std::cout << "FluidSystem::precipitateMolarDensity("<<phaseIdx<<") = " << FluidSystem::precipitateMolarDensity(phaseIdx) << "\n";
//             for (int phaseIdx=1; phaseIdx<numSPhases; ++phaseIdx)
            return FluidSystem::precipitateMolarDensity(phaseIdx);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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
    { return this->fluidState_.moleFraction(phaseIdx, compIdx)
                  /(this->fluidState_.moleFraction(phaseIdx, phaseIdx)
                  * FluidSystem::molarMass(phaseIdx));}

   /*!
    * Circumvents the inheritance architecture of the nonisothermal model
    */
    static Scalar callProtectedTemperature(const ElementSolutionVector& elemSol,
                                        const Problem& problem,
                                        const Element& element,
                                        const SubControlVolume& scv)
    {
//          return Implementation::temperature_(priVars, problem,element, fvGeometry, scvIdx);
        return BaseType::temperature(elemSol, problem, element, scv);
    }

    /*!
    * Circumvents the inheritance architecture of the ninisothermal model
    */
    void callProtectedUpdateEnergy(const ElementSolutionVector& elemSol,
                                        const Problem& problem,
                                        const Element& element,
                                        const SubControlVolume& scv)
    {
        asImp_().updateEnergy_(elemSol, problem, element, scv);
    };

protected:
    friend class OnePNCVolumeVariables<TypeTag>;
    static Scalar temperature_(const ElementSolutionVector& elemSol,
                                        const Problem& problem,
                                        const Element& element,
                                        const SubControlVolume& scv)
    {
        return problem.temperatureAtPos(scv);
    }

    template<class ParameterCache>
    static Scalar enthalpy_(const FluidState& fluidState,
                            const ParameterCache& paramCache,
                            int phaseIdx)
    {
        return 0;
    }

   /*!
    * \brief Update all quantities for a given control volume.
    *
    * \param priVars The solution primary variables
    * \param problem The problem
    * \param element The element
    * \param fvGeometry Evaluate function with solution of current or previous time step
    * \param scvIdx The local index of the SCV (sub-control volume)
    * \param isOldSol Evaluate function with solution of current or previous time step
    */
    void updateEnergy_(const ElementSolutionVector& elemSol,
                       const Problem& problem,
                       const Element& element,
                       const SubControlVolume& scv)
    { };

    Scalar precipitateVolumeFraction_[numSPhases];
    Scalar permeabilityFactor_;
    Scalar initialPorosity_;
    Scalar initialPermeability_;
    Scalar minimumPorosity_;
    Scalar maximumPorosity_;
    Scalar sumPrecipitates_;


private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace

#endif
