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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase n-component mineralisation model.
 */

#ifndef DUMUX_1PNCMIN_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_1PNCMIN_LOCAL_RESIDUAL_BASE_HH

#include "properties.hh"
#include <dumux/porousmediumflow/compositional/localresidual.hh>

namespace Dumux
{
/*!
 * \ingroup OnePNCMinModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase n-component mineralization fully implicit model.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the one-phase n-component flow.
 */
template<class TypeTag>
class OnePNCMinLocalResidual: public CompositionalLocalResidual<TypeTag>
{
protected:
//     typedef OnePNCLocalResidual<TypeTag> ParentType;
    using ParentType = CompositionalLocalResidual<TypeTag>;
    using ThisType = OnePNCMinLocalResidual<TypeTag>;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    enum
    {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        pressureIdx = Indices::pressureIdx,
        firstMoleFracIdx = Indices::firstMoleFracIdx,

        phaseIdx = Indices::phaseIdx,
        sPhaseIdx = 1, // don't use the sPhaseIdx of the fluidsystem

        conti0EqIdx = Indices::conti0EqIdx,
    };


public:
    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume).
     * In contrast to the 1pnc model, here, the storage of solid phases is included too.
     *
     *  \param scv the SCV (sub-control-volume)
     *  \param volVars The volume variables of the right time step
     */
    PrimaryVariables computeStorage(const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        //call parenttype function
        auto storage = ParentType::computeStorage(scv, volVars);

        // Compute storage term of all solid (precipitated) phases (excluding the non-reactive matrix)
         for (int Idx = sPhaseIdx; Idx < numPhases + numSPhases; ++Idx)
        {
            auto eqIdx = conti0EqIdx + numComponents-numPhases + Idx;
            storage[eqIdx] += volVars.precipitateVolumeFraction(Idx)*volVars.molarDensity(Idx);
        }

        return storage;
    }
};

} // end namespace Dumux

#endif
