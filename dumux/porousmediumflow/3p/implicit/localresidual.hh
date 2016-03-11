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
 *        using the three-phase fully implicit model.
 */
#ifndef DUMUX_3P_LOCAL_RESIDUAL_HH
#define DUMUX_3P_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase fully implicit model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for three-phase flow.
 */
template<class TypeTag>
class ThreePLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        conti0EqIdx = Indices::conti0EqIdx //!< Index of the mass conservation equation for the water component
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

public:
    ThreePLocalResidual() : ParentType()
    {
        massWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }
    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub-control volume divided by the volume).
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    PrimaryVariables computeStorage(const SubControlVolume &scv,
                                    const VolumeVariables &volVars) const
    {
        PrimaryVariables storage(0.0);

        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            auto eqIdx = conti0EqIdx + phaseIdx;
            storage[eqIdx] +=
                volVars.porosity()
                * volVars.saturation(phaseIdx)
                * volVars.density(phaseIdx);
        }

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param fIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    PrimaryVariables computeFlux(const SubControlVolumeFace& scvFace) const
    {
        auto& fluxVars = this->model_().fluxVars(scvFace);

        // get upwind weights into local scope
        auto massWeight = massWeight_;
        PrimaryVariables flux(0.0);

        // advective fluxes of all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            auto upwindRule = [massWeight, phaseIdx](const VolumeVariables& up, const VolumeVariables& dn)
            {
                return (massWeight)*up.mobility(phaseIdx)*up.density(phaseIdx)
                      +(1-massWeight)*dn.mobility(phaseIdx)*dn.density(phaseIdx);
            };

            // add advective flux of current phase
            auto eqIdx = conti0EqIdx + phaseIdx;
            flux[eqIdx] += fluxVars.advection().flux(phaseIdx, upwindRule);
        }

        return flux;
    }

protected:
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }

    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }
private:
    Scalar massWeight_;
};

} // end namespace

#endif
