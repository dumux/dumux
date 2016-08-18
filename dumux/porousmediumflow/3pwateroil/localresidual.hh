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
#ifndef DUMUX_3P2CNI_LOCAL_RESIDUAL_HH
#define DUMUX_3P2CNI_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePWaterOilModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase three-component fully implicit model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the 3P2C flow.
 */
template<class TypeTag>
class ThreePWaterOilLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;


    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        conti0EqIdx = Indices::conti0EqIdx,//!< Index of the mass conservation equation for the water component
        conti1EqIdx = Indices::conti1EqIdx,//!< Index of the mass conservation equation for the contaminant component
        energyEqIdx = Indices::energyEqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
public:

    /*!
     * \brief Evaluate the storage term of the current solution in a
     *        single phase.
     *
     * \param element The element
     * \param phaseIdx The index of the fluid phase
     */
    void evalPhaseStorage(const Element &element, const int phaseIdx)
    {
        FVElementGeometry fvGeometry;
        fvGeometry.update(this->gridView_(), element);
        ElementBoundaryTypes bcTypes;
        bcTypes.update(this->problem_(), element, fvGeometry);
        ElementVolumeVariables elemVolVars;
        elemVolVars.update(this->problem_(), element, fvGeometry, false);

        this->storageTerm_.resize(fvGeometry.numScv);
        this->storageTerm_ = 0;

        this->elemPtr_ = &element;
        this->fvElemGeomPtr_ = &fvGeometry;
        this->bcTypesPtr_ = &bcTypes;
        this->prevVolVarsPtr_ = 0;
        this->curVolVarsPtr_ = &elemVolVars;
        evalPhaseStorage_(phaseIdx);
    }

    /*!
     * \brief Evaluate the amount all conservation quantities
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
                    * volVars.fluidState().moleFraction(phaseIdx, compIdx);
            }
        }

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.

        // compute the energy storage
        storage[energyEqIdx] = volVars.porosity()
            *(
                volVars.density(wPhaseIdx)
                *volVars.internalEnergy(wPhaseIdx)
                *volVars.saturation(wPhaseIdx)
                +
                volVars.density(nPhaseIdx)
                *volVars.internalEnergy(nPhaseIdx)
                *volVars.saturation(nPhaseIdx)
                +
                volVars.density(gPhaseIdx)
                *volVars.internalEnergy(gPhaseIdx)
                *volVars.saturation(gPhaseIdx)
            )
            + volVars.temperature()*volVars.heatCapacity();
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param faceIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, const bool onBoundary=false) const
    {
        FluxVariables fluxVars;
        fluxVars.update(this->problem_(),
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
                       * up.fluidState().molarDensity(phaseIdx)
                       * up.fluidState().moleFraction(phaseIdx, compIdx)
                       +
                       (1.0 - massUpwindWeight)
                       * dn.fluidState().molarDensity(phaseIdx)
                       * dn.fluidState().moleFraction(phaseIdx, compIdx));
            }
        }

        // advective heat flux in all phases
        flux[energyEqIdx] = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // vertex data of the upstream and the downstream vertices
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

            flux[energyEqIdx] +=
                fluxVars.volumeFlux(phaseIdx) * (
                                              massUpwindWeight * // upstream vertex
                                              (  up.density(phaseIdx) *
                                                 up.enthalpy(phaseIdx))
                                              +
                                              (1-massUpwindWeight) * // downstream vertex
                                              (  dn.density(phaseIdx) *
                                                 dn.enthalpy(phaseIdx)) );
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
        Scalar tmp;
        tmp = - fluxVars.porousDiffCoeff()[wPhaseIdx] * fluxVars.molarDensity(wPhaseIdx);
        tmp *= (fluxVars.moleFractionCompNGrad(wPhaseIdx) * fluxVars.face().normal);
        Scalar jNW = tmp;

        Scalar jWW = -jNW;

        tmp = - fluxVars.porousDiffCoeff()[gPhaseIdx] * fluxVars.molarDensity(gPhaseIdx);
        tmp *= (fluxVars.moleFractionCompWGrad(gPhaseIdx) * fluxVars.face().normal);
        Scalar jWG = tmp;

        Scalar jNG = -jWG;

        tmp = - fluxVars.porousDiffCoeff()[nPhaseIdx] * fluxVars.molarDensity(nPhaseIdx);
        tmp *= (fluxVars.moleFractionCompWGrad(nPhaseIdx) * fluxVars.face().normal);
        Scalar jWN = tmp;

        Scalar jNN = -jWN;

        flux[conti0EqIdx] += jWW+jWG+jWN;
        flux[conti1EqIdx] += jNW+jNG+jNN;
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
    void evalPhaseStorage_(const int phaseIdx)
    {
        if(!useMoles) //mass-fraction formulation
        {
            // evaluate the storage terms of a single phase
            for (int i=0; i < this->fvGeometry_().numScv; i++) {
                PrimaryVariables &storage = this->storageTerm_[i];
                const ElementVolumeVariables &elemVolVars = this->curVolVars_();
                const VolumeVariables &volVars = elemVolVars[i];

                // compute storage term of all components within all phases
                storage = 0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    int eqIdx = (compIdx == wCompIdx) ? conti0EqIdx : conti1EqIdx;
                    storage[eqIdx] += volVars.density(phaseIdx)
                        * volVars.saturation(phaseIdx)
                        * volVars.fluidState().massFraction(phaseIdx, compIdx);
                }

                storage *= volVars.porosity();
                storage *= this->fvGeometry_().subContVol[i].volume;
            }
        }
        else //mole-fraction formulation
        {
            // evaluate the storage terms of a single phase
            for (int i=0; i < this->fvGeometry_().numScv; i++) {
                PrimaryVariables &storage = this->storageTerm_[i];
                const ElementVolumeVariables &elemVolVars = this->curVolVars_();
                const VolumeVariables &volVars = elemVolVars[i];

                // compute storage term of all components within all phases
                storage = 0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    int eqIdx = (compIdx == wCompIdx) ? conti0EqIdx : conti1EqIdx;
                    storage[eqIdx] += volVars.molarDensity(phaseIdx)
                        * volVars.saturation(phaseIdx)
                        * volVars.fluidState().moleFraction(phaseIdx, compIdx);
                }

                storage *= volVars.porosity();
                storage *= this->fvGeometry_().subContVol[i].volume;
            }
        }
    }
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
