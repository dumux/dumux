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
 *        using the three-phase three-component fully implicit model.
 */
#ifndef DUMUX_3P3C_LOCAL_RESIDUAL_HH
#define DUMUX_3P3C_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase three-component fully implicit model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the 3P3C flow.
 */
template<class TypeTag>
class ThreePThreeCLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        conti0EqIdx = Indices::conti0EqIdx,//!< Index of the mass conservation equation for the water component
        conti1EqIdx = Indices::conti1EqIdx,//!< Index of the mass conservation equation for the contaminant component
        conti2EqIdx = Indices::conti2EqIdx,//!< Index of the mass conservation equation for the gas component

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
        gCompIdx = Indices::gCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

public:

    ThreePThreeCLocalResidual() : ParentType()
    {
        massWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }
    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    PrimaryVariables computeStorage(const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        PrimaryVariables storage(0.0);

        // compute storage term of all components within all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                auto eqIdx = conti0EqIdx + compIdx;
                storage[eqIdx] += volVars.porosity()
                                  * volVars.saturation(phaseIdx)
                                  * volVars.molarDensity(phaseIdx)
                                  * volVars.moleFraction(phaseIdx, compIdx);
            }
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
    PrimaryVariables computeFlux(const SubControlVolumeFace& scvFace)
    {
        auto& fluxVars = this->model_().fluxVars_(scvFace);

        // get upwind weights into local scope
        auto massWeight = massWeight_;
        PrimaryVariables flux(0.0);

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                auto upwindRule = [massWeight, phaseIdx, compIdx](const VolumeVariables& up, const VolumeVariables& dn)
                {
                    return ( massWeight )*up.molarDensity(phaseIdx)
                                         *up.moleFraction(phaseIdx, compIdx)
                                         *up.mobility(phaseIdx)
                          +(1-massWeight)*dn.molarDensity(phaseIdx)
                                         *dn.moleFraction(phaseIdx, compIdx)
                                         *dn.mobility(phaseIdx);
                };

                // get equation index
                auto eqIdx = conti0EqIdx + compIdx;
                flux[eqIdx] += fluxVars.advection().flux(phaseIdx, upwindRule);
            }
        }

        // diffusive fluxes
        Scalar jGW = fluxVars.molecularDiffusion(wPhaseIdx, gCompIdx).flux();
        Scalar jNW = fluxVars.molecularDiffusion(wPhaseIdx, nCompIdx).flux();
        Scalar jWW = -(jGW+jNW);

        Scalar jWG = fluxVars.molecularDiffusion(gPhaseIdx, wCompIdx).flux();
        Scalar jNG = fluxVars.molecularDiffusion(gPhaseIdx, nCompIdx).flux();
        Scalar jGG = -(jWG+jNG);

        // Scalar jWN = fluxVars.molecularDiffusion().flux(nPhaseIdx, wCompIdx);
        // Scalar jGN = fluxVars.molecularDiffusion().flux(nPhaseIdx, gCompIdx);
        // Scalar jNN = -(jGN+jWN);

        // At the moment we do not consider diffusion in the NAPL phase
        Scalar jWN = 0;
        Scalar jGN = 0;
        Scalar jNN = 0;

        flux[conti0EqIdx] += jWW+jWG+jWN;
        flux[conti1EqIdx] += jNW+jNG+jNN;
        flux[conti2EqIdx] += jGW+jGG+jGN;

        return flux;
    }

protected:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

private:
    Scalar massWeight_;
};

} // end namespace

#endif
