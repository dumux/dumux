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
 *
 * \brief Element-wise calculation of the local residual for problems
 *        using the two-phase n-component mineralisation box model.
 */

#ifndef DUMUX_2PNCMIN_LOCAL_RESIDUAL_HH
#define DUMUX_2PNCMIN_LOCAL_RESIDUAL_HH

#include "properties.hh"
#include <dumux/porousmediumflow/compositional/localresidual.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPNCMinModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the local residual for problems
 *        using the two-phase n-component mineralization fully implicit model.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the two-phase n-component flow.
 */
template<class TypeTag>
class TwoPNCMinLocalResidual : public CompositionalLocalResidual<TypeTag>
{
    using ParentType = CompositionalLocalResidual<TypeTag>;
    using ThisType = TwoPNCMinLocalResidual<TypeTag>;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        conti0EqIdx = Indices::conti0EqIdx
    };

public:
    using ParentType::ParentType;
    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume).
     * In contrast to the 2pnc model, here, the storage of solid phases is included too.
     *
     *  \param scv the SCV (sub-control-volume)
     *  \param volVars The volume variables of the right time step
     */
    PrimaryVariables computeStorage(const Problem& problem,
                                    const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        // call parenttype function
        auto storage = ParentType::computeStorage(problem, scv, volVars);

        // Compute storage term of all solid (precipitated) phases (excluding the non-reactive matrix)
        for (int phaseIdx = numPhases; phaseIdx < numPhases + numSPhases; ++phaseIdx)
        {
            auto eqIdx = conti0EqIdx + numComponents-numPhases + phaseIdx;
            storage[eqIdx] += volVars.precipitateVolumeFraction(phaseIdx)*volVars.molarDensity(phaseIdx);
        }

        return storage;
    }
};

} // end namespace Dumux

#endif
