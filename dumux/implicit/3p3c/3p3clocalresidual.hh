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

#include "3p3cproperties.hh"

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
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
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
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

public:
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
    void computeStorage(PrimaryVariables &storage, const int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars =
            usePrevSol
            ? this->prevVolVars_()
            : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute storage term of all components within all phases
        storage = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                storage[conti0EqIdx + compIdx] +=
                    volVars.porosity()
                    * volVars.saturation(phaseIdx)
                    * volVars.molarDensity(phaseIdx)
                    * volVars.moleFraction(phaseIdx, compIdx);
            }
        }
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
    void computeFlux(PrimaryVariables &flux, const int fIdx, const bool onBoundary=false) const
    {
        FluxVariables fluxVars(this->problem_(),
                           this->element_(),
                           this->fvGeometry_(),
                           fIdx,
                           this->curVolVars_(),
                           onBoundary);

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */

    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        Scalar massUpwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);

        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // add advective flux of current component in current
                // phase
                // if alpha > 0 und alpha < 1 then both upstream and downstream
                // nodes need their contribution
                // if alpha == 1 (which is mostly the case) then, the downstream
                // node is not evaluated
                int eqIdx = conti0EqIdx + compIdx;
                flux[eqIdx] += fluxVars.volumeFlux(phaseIdx)
                    * (massUpwindWeight
                       * up.molarDensity(phaseIdx)
                       * up.moleFraction(phaseIdx, compIdx)
                       +
                       (1.0 - massUpwindWeight)
                       * dn.molarDensity(phaseIdx)
                       * dn.moleFraction(phaseIdx, compIdx));
            }
        }
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */

    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // TODO: reference!?  Dune::FieldMatrix<Scalar, numPhases, numComponents> averagedPorousDiffCoeffMatrix = fluxVars.porousDiffCoeff();
        // add diffusive flux of gas component in liquid phase
        Scalar tmp = - fluxVars.porousDiffCoeff()[wPhaseIdx][gCompIdx] * fluxVars.molarDensity(wPhaseIdx);
        tmp *= (fluxVars.moleFractionCompGGrad(wPhaseIdx) * fluxVars.face().normal);
        Scalar jGW = tmp;

        tmp = - fluxVars.porousDiffCoeff()[wPhaseIdx][nCompIdx] * fluxVars.molarDensity(wPhaseIdx);
        tmp *= (fluxVars.moleFractionCompNGrad(wPhaseIdx) * fluxVars.face().normal);
        Scalar jNW = tmp;

        Scalar jWW = -(jGW+jNW);

        tmp = - fluxVars.porousDiffCoeff()[gPhaseIdx][wCompIdx] * fluxVars.molarDensity(gPhaseIdx);
        tmp *= (fluxVars.moleFractionCompWGrad(gPhaseIdx) * fluxVars.face().normal);
        Scalar jWG = tmp;

        tmp = - fluxVars.porousDiffCoeff()[gPhaseIdx][nCompIdx] * fluxVars.molarDensity(gPhaseIdx);
        tmp *= (fluxVars.moleFractionCompNGrad(gPhaseIdx) * fluxVars.face().normal);
        Scalar jNG = tmp;

        Scalar jGG = -(jWG+jNG);

        tmp = - fluxVars.porousDiffCoeff()[nPhaseIdx][wCompIdx] * fluxVars.molarDensity(nPhaseIdx);
        tmp *= (fluxVars.moleFractionCompWGrad(nPhaseIdx) * fluxVars.face().normal);
        Scalar jWN = tmp;

        tmp = - fluxVars.porousDiffCoeff()[nPhaseIdx][gCompIdx] * fluxVars.molarDensity(nPhaseIdx);
        tmp *= (fluxVars.moleFractionCompGGrad(nPhaseIdx) * fluxVars.face().normal);
        Scalar jGN = tmp;

        Scalar jNN = -(jGN+jWN);

        flux[conti0EqIdx] += jWW+jWG+jWN;
        flux[conti1EqIdx] += jNW+jNG+jNN;
        flux[conti2EqIdx] += jGW+jGG+jGN;
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source The source/sink in the SCV for each component
     * \param scvIdx The index of the SCV
     */
    void computeSource(PrimaryVariables &source, const int scvIdx)
    {
        this->problem_().solDependentSource(source,
                                     this->element_(),
                                     this->fvGeometry_(),
                                     scvIdx,
                                     this->curVolVars_());
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
};

} // end namespace

#endif
