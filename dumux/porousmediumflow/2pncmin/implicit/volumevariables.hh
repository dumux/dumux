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
    // base type is used for energy related quantites
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

        useSalinity = GET_PROP_VALUE(TypeTag, useSalinity)
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CoordScalar = typename Grid::ctype;
    using Miscible2pNCComposition = Dumux::Miscible2pNCComposition<Scalar, FluidSystem>;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FluidSystem>;

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
        // Update parent type (also completes the fluid state)
        ParentType::update(elemSol, problem, element, scv);

        Implementation::completeFluidState(elemSol, problem, element, scv, fluidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        auto&& priVars = isBox ? elemSol[scv.indexInElement()] : elemSol[0];

        sumPrecipitates_ = 0.0;
        for(int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
        {
           precipitateVolumeFraction_[sPhaseIdx] = priVars[numComponents + sPhaseIdx];
           sumPrecipitates_+= precipitateVolumeFraction_[sPhaseIdx];
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

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // end namespace

#endif
