// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-     by Holger Class                                 *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 *        using the non-isothermal three-phase three-component box model.
 *
 */
#ifndef DUMUX_NEW_3P3CNI_LOCAL_RESIDUAL_HH
#define DUMUX_NEW_3P3CNI_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/3p3c/3p3clocalresidual.hh>


#include "3p3cnivolumevariables.hh"
#include "3p3cnifluxvariables.hh"

#include "3p3cniproperties.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCNIModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase three-component box model.
 */
template<class TypeTag>
class ThreePThreeCNILocalResidual : public ThreePThreeCLocalResidual<TypeTag>
{
    typedef ThreePThreeCLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, ThreePThreeCNIIndices) Indices;

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


    static constexpr Scalar massUpwindWeight = GET_PROP_VALUE(TypeTag, MassUpwindWeight);

public:
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
        Scalar wdens = vertDat.density(wPhaseIdx);
        Scalar ndens = vertDat.density(nPhaseIdx);
        Scalar gdens = vertDat.density(gPhaseIdx);
        Scalar wInEnerg = vertDat.internalEnergy(wPhaseIdx);
        Scalar nInEnerg = vertDat.internalEnergy(nPhaseIdx);
        Scalar gInEnerg = vertDat.internalEnergy(gPhaseIdx);
        Scalar wsat = vertDat.saturation(wPhaseIdx);
        Scalar nsat = vertDat.saturation(nPhaseIdx);
        Scalar gsat = vertDat.saturation(gPhaseIdx);
        Scalar temp = vertDat.temperature();
        Scalar heatCap = vertDat.heatCapacity();
        Scalar poro = vertDat.porosity();

        result[energyEqIdx] = temp*heatCap
            + poro * (gdens*gInEnerg*gsat
                      + wdens*wInEnerg*wsat
                      + ndens*nInEnerg*nsat);
        /*
          vertDat.porosity()*(vertDat.density(wPhaseIdx) *
          vertDat.internalEnergy(wPhaseIdx) *
          vertDat.saturation(wPhaseIdx)
          +
          vertDat.density(nPhaseIdx) *
          vertDat.internalEnergy(nPhaseIdx) *
          vertDat.saturation(nPhaseIdx)
          +
          vertDat.density(gPhaseIdx) *
          vertDat.internalEnergy(gPhaseIdx) *
          vertDat.saturation(gPhaseIdx))
          +
          vertDat.temperature()*vertDat.heatCapacity();
        */
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
                                              massUpwindWeight * // upstream vertex
                                              (  up.density(phase) *
                                                 up.mobility(phase) *
                                                 up.enthalpy(phase))
                                              +
                                              (1-massUpwindWeight) * // downstream vertex
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
};

}

#endif
