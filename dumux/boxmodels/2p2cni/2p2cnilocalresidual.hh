// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal two-phase two-component box model.
 *
 */
#ifndef DUMUX_NEW_2P2CNI_LOCAL_RESIDUAL_HH
#define DUMUX_NEW_2P2CNI_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/2p2c/2p2clocalresidual.hh>

#include <dumux/boxmodels/2p2cni/2p2cnivolumevariables.hh>
#include <dumux/boxmodels/2p2cni/2p2cnifluxvariables.hh>
#include <dumux/boxmodels/2p2cni/2p2cniproperties.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCNIModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component box model.
 */
template<class TypeTag>
class TwoPTwoCNILocalResidual : public TwoPTwoCLocalResidual<TypeTag>
{
    typedef TwoPTwoCLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryVariables) BoundaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        energyEqIdx = Indices::energyEqIdx,
        temperatureIdx = Indices::temperatureIdx,
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx
    };

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    TwoPTwoCNILocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param result The storage of the conservation quantitiy (mass or energy) within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &result, int scvIdx, bool usePrevSol) const
    {
        // compute the storage term for phase mass
        ParentType::computeStorage(result, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemDat = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &vertDat = elemDat[scvIdx];

        // compute the energy storage
        result[energyEqIdx] =
            vertDat.porosity()*(vertDat.density(lPhaseIdx) *
                                vertDat.internalEnergy(lPhaseIdx) *
                                //vertDat.enthalpy(lPhaseIdx) *
                                vertDat.saturation(lPhaseIdx)
                                +
                                vertDat.density(gPhaseIdx) *
                                vertDat.internalEnergy(gPhaseIdx) *
                                //vertDat.enthalpy(gPhaseIdx) *
                                vertDat.saturation(gPhaseIdx))
            +
            vertDat.temperature()*vertDat.heatCapacity();
    }

    /*!
     * \brief Evaluates the advective mass flux and the heat flux
     * over a face of a subcontrol volume and writes the result in
     * the flux vector.
     *
     * \param flux The advective flux over the SCV (sub-control-volume) face for each component
     * \param fluxData The flux variables at the current SCV face
     *
     * This method is called by compute flux (base class)
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxData) const
    {
        // advective mass flux
        ParentType::computeAdvectiveFlux(flux, fluxData);

        // advective heat flux in all phases
        flux[energyEqIdx] = 0;
        for (int phase = 0; phase < numPhases; ++phase) {
            // vertex data of the upstream and the downstream vertices
            const VolumeVariables &up = this->curVolVars_(fluxData.upstreamIdx(phase));
            const VolumeVariables &dn = this->curVolVars_(fluxData.downstreamIdx(phase));

            flux[energyEqIdx] +=
                fluxData.KmvpNormal(phase) * (
                    massUpwindWeight_ * // upstream vertex
                    (  up.density(phase) *
                       up.mobility(phase) *
                       up.enthalpy(phase))
                    +
                    (1-massUpwindWeight_) * // downstream vertex
                    (  dn.density(phase) *
                       dn.mobility(phase) *
                       dn.enthalpy(phase)) );
        }
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the SCV (sub-control-volume) face for each conservation quantity (mass, energy)
     * \param fluxData The flux variables at the current SCV face
     *
     * This method is called by compute flux (base class)
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxData) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxData);

        // diffusive heat flux
        flux[temperatureIdx] +=
            fluxData.normalMatrixHeatFlux();
    }

    /*!
     * \brief Compute the fluxes at outflow boundaries. This does essentially the same
     *        as computeFluxes, but the fluxes are evaluated at the integration point
     *        of the boundary face. Therefore, still some variables are evaluated
     *        at the vertex (usually the ones which are upwinded).
     *
     * \param flux A temporary vector, where the outflow boundary fluxes are written into
     * \param boundaryVars The boundary variables object
     * \param scvIdx The index of the SCV containing the outflow boundary face
     * \param boundaryFaceIdx The index of the boundary face
     */
    void computeOutflowValues(PrimaryVariables &flux,
            const BoundaryVariables &boundaryVars,
            const int scvIdx,
            const int boundaryFaceIdx)
    {
        ParentType::computeOutflowValues(flux, boundaryVars, scvIdx, boundaryFaceIdx);

        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
        const VolumeVariables &vertVars = this->curVolVars_()[scvIdx];

        if (bcTypes.isOutflow(energyEqIdx))
        {
            // advective heat flux in all phases
            flux[energyEqIdx] = 0;
            for (int phase = 0; phase < numPhases; ++phase) {
                flux[energyEqIdx] +=
                    boundaryVars.KmvpNormal(phase) *
                            vertVars.density(phase) *
                            vertVars.mobility(phase) *
                            vertVars.enthalpy(phase);
            }

            // add conductive heat flux in all phases
            flux[energyEqIdx] +=
                boundaryVars.normalMatrixHeatFlux();
        }
    }

private:
    Scalar massUpwindWeight_;

};

}

#endif
