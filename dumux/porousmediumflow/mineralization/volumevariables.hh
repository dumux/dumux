// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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

namespace Dumux {

/*!
 * \ingroup MineralizationModel
 * \brief Contains the quantities which are are constant within a sub-control volume
 *        of the finite volume grid in an m-phase, n-component, mineralization model.
 */
template <class Traits, class NonMineralizationVolVars>
class MineralizationVolumeVariables : public NonMineralizationVolVars
{
    using ParentType = NonMineralizationVolVars;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;
    using SolidState = typename Traits::SolidState;

public:
    using SolidSystem = typename Traits::SolidSystem;

    /*!
     * \brief Returns the volume fraction of the precipitate (solid phase)
     *        for the given phaseIdx.
     *
     * \param sCompIdx The index of the solid component
     */
    Scalar solidVolumeFraction(int sCompIdx) const
    { return this->solidState_.volumeFraction(sCompIdx); }

    /*!
     * \brief Returns the density of the phase for all fluid and solid phases.
     *
     * \param sCompIdx The index of the solid component
     */
    Scalar solidComponentDensity(int sCompIdx) const
    {
        return SolidSystem::density(this->solidState_, sCompIdx);
    }

    /*!
     * \brief Returns the density of the phase for all fluid and solid phases.
     *
     * \param sCompIdx The index of the solid component
     */
    Scalar solidComponentMolarDensity(int sCompIdx) const
    {
        return SolidSystem::molarDensity(this->solidState_, sCompIdx);
    }

};
} // end namespace Dumux

#endif
