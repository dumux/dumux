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
 * \brief Element-wise calculation of the residual for problems
 *        using the non-isothermal three-phase three-component fully implicit model.
 *
 */
#ifndef DUMUX_NEW_3P3CNI_LOCAL_RESIDUAL_HH
#define DUMUX_NEW_3P3CNI_LOCAL_RESIDUAL_HH

#include "3p3cniproperties.hh"
#include <dumux/implicit/3p3c/3p3clocalresidual.hh>

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCNIModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the residual for problems
 *        using the non-isothermal three-phase three-component fully implicit model.
 */
template<class TypeTag>
class ThreePThreeCNILocalResidual : public ThreePThreeCLocalResidual<TypeTag>
{
    typedef ThreePThreeCLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        energyEqIdx = Indices::energyEqIdx,
        temperatureIdx = Indices::temperatureIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx
    };


    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    
public:
    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The storage of the conservation quantitiy (mass or energy) within the sub-control volume
     *  \param scvIdx The sub-control-volume index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, bool usePrevSol) const
    {
        // compute the storage term for phase mass
        ParentType::computeStorage(storage, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

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
     * \brief Evaluates the advective mass flux and the heat flux
     *        over a face of a subcontrol volume and writes the result in
     *        the flux vector.
     *
     * \param flux The advective flux over the SCV (sub-control-volume) face for each component
     * \param fluxVars The flux variables at the current SCV face
     *
     * This method is called by compute flux (base class).
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // advective mass flux
        ParentType::computeAdvectiveFlux(flux, fluxVars);

        static const Scalar massUpwindWeight = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);

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
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each conservation quantity (mass, energy)
     * \param fluxVars The flux variables at the current sub-control-volume face
     *
     * This method is called by compute flux (base class).
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxVars);

        // diffusive heat flux
        flux[temperatureIdx] +=
            fluxVars.normalMatrixHeatFlux();
    }
};

}

#endif
