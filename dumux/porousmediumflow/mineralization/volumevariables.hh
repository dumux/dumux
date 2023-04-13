// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
