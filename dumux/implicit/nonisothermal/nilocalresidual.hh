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
 *        using the non-isothermal two-phase two-component fully implicit model.
 *
 */
#ifndef DUMUX_NEW_NI_LOCAL_RESIDUAL_HH
#define DUMUX_NEW_NI_LOCAL_RESIDUAL_HH

#include "niproperties.hh"

namespace Dumux
{
/*!
 * \ingroup NIModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component fully implicit model.
 */
template<class TypeTag>
class NILocalResidual : public GET_PROP_TYPE(TypeTag, IsothermalLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, IsothermalLocalResidual) ParentType;


    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        energyEqIdx = Indices::energyEqIdx,
        temperatureIdx = Indices::temperatureIdx,

    };

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    NILocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The storage of the conservation quantity (mass or energy) within the sub-control volume
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
        // using the implicit Euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];
        storage[energyEqIdx] = 0;

        // add the contribution from the fluid phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            storage[energyEqIdx] +=
                volVars.porosity()*(volVars.fluidState().density(phaseIdx) *
                                    volVars.fluidState().internalEnergy(phaseIdx) *
                                    volVars.fluidState().saturation(phaseIdx));

        }

        // add the contribution from the solid phase
        storage[energyEqIdx] += volVars.temperature()
                               *volVars.solidHeatCapacity()
                               *volVars.solidDensity()
                               *(1 - volVars.porosity());
    }

    /*!
     * \brief Evaluates the advective mass flux and the heat flux
     *        over a face of a sub-control volume and writes the result into
     *        the flux vector.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current sub-control-volume face
     *
     * This method is called by compute flux (base class).
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // advective mass flux
        ParentType::computeAdvectiveFlux(flux, fluxVars);

        // advective heat flux in all phases
        flux[energyEqIdx] = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // vertex data of the upstream and the downstream vertices
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(phaseIdx));

            flux[energyEqIdx] +=
                fluxVars.volumeFlux(phaseIdx) * (massUpwindWeight_ * // upstream vertex
                                                 (  up.fluidState().density(phaseIdx) *
                                                    up.fluidState().enthalpy(phaseIdx))
                                                 +
                                                 (1-massUpwindWeight_) * // downstream vertex
                                                 (  dn.fluidState().density(phaseIdx) *
                                                    dn.fluidState().enthalpy(phaseIdx)) );

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

private:
    Scalar massUpwindWeight_;



};

}

#endif
