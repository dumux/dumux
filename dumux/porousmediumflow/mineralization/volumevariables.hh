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
 * \ingroup MineralizationModel
 * \brief Contains the quantities which are constant within a sub-control volume
 *        of the finite volume grid in the two-phase, n-component mineralization model.
 */
#ifndef DUMUX_MINERALIZATION_VOLUME_VARIABLES_HH
#define DUMUX_MINERALIZATION_VOLUME_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{

/*!
 * \ingroup MineralizationModel
 * \brief Contains the quantities which are are constant within a sub-control volume
 *        of the finite volume grid in an m-phase, n-component, mineralization model.
 */
template <class TypeTag>
class MineralizationVolumeVariables : public GET_PROP_TYPE(TypeTag, NonMineralizationVolumeVariables)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, NonMineralizationVolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Element = typename GridView::template Codim<0>::Entity;

    enum
    {
        dimWorld=GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases =  GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
    };

public:
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    //! updates all required quantities inside the given scv
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        // Update parent type (also completes the fluid state)
        ParentType::update(elemSol, problem, element, scv);

        // calculate the remaining quantities
        auto&& priVars = elemSol[scv.indexInElement()];

        sumPrecipitates_ = 0.0;
        for(int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
        {
           precipitateVolumeFraction_[sPhaseIdx] = priVars[numComponents + sPhaseIdx];
           sumPrecipitates_+= precipitateVolumeFraction_[sPhaseIdx];
        }
    }

    /*!
     * \brief Returns the volume fraction of the precipitate (solid phase)
     *        for the given phaseIdx
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
     * \f$\mathrm{molality}
     *  = \frac{n_\mathrm{component}}{m_\mathrm{solvent}}
     *  = \frac{n_\mathrm{component}}{n_\mathrm{solvent}*M_\mathrm{solvent}}\f$
     *
     * \note compIdx of the main component (solvent) in the
     *       phase is equal to the phaseIdx
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
};
} // end namespace Dumux

#endif
