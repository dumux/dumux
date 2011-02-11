// $Id: 2pnilocalresidual.hh 3840 2010-07-15 10:14:15Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 *        using the non-isothermal two-phase box model.
 *
 */
#ifndef DUMUX_NEW_2PNI_LOCAL_RESIDUAL_HH
#define DUMUX_NEW_2PNI_LOCAL_RESIDUAL_HH

#include "2pniproperties.hh"

#include <dumux/boxmodels/2p/2plocalresidual.hh>


#include <dumux/boxmodels/2pni/2pnivolumevariables.hh>
#include <dumux/boxmodels/2pni/2pnifluxvariables.hh>


namespace Dumux
{
/*!
 * \ingroup TwoPNIBoxModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal two-phase box model.
 */

template<class TypeTag>
class TwoPNILocalResidual : public TwoPLocalResidual<TypeTag>
{
    typedef TwoPNILocalResidual<TypeTag> ThisType;
    typedef TwoPLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

    static const Scalar mobilityUpwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass and energy storage) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param result The phase mass within the sub-control volume
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
        const ElementVolumeVariables &vertDatArray = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &vertDat = vertDatArray[scvIdx];

        // compute the energy storage
        result[temperatureIdx] =
            vertDat.porosity()*(vertDat.density(wPhaseIdx) *
                                vertDat.internalEnergy(wPhaseIdx) *
                                vertDat.saturation(wPhaseIdx)
                                +
                                vertDat.density(nPhaseIdx) *
                                vertDat.internalEnergy(nPhaseIdx) *
                                vertDat.saturation(nPhaseIdx))
          + vertDat.temperature()*vertDat.heatCapacity();
    }

    /*!
     * \brief Evaluates the advective mass flux and the heat flux
     * over a face of a subcontrol volume and writes the result in
     * the flux vector.
     *
     * \param flux The advective flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     *
     * This method is called by compute flux (base class)
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // advective mass flux
        ParentType::computeAdvectiveFlux(flux, fluxVars);

        // advective heat flux in all phases
        flux[energyEqIdx] = 0;
        Vector tmpVec;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // calculate the flux in the normal direction of the
            // current sub control volume face
            fluxVars.intrinsicPermeability().mv(fluxVars.potentialGrad(phaseIdx),
                                                tmpVec);
            Scalar normalFlux = - (tmpVec*fluxVars.face().normal);

            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(normalFlux));
            const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(normalFlux));

            // add advective energy flux in current phase
            flux[energyEqIdx] +=
                normalFlux
                *
                ((    mobilityUpwindAlpha)*
                 up.density(phaseIdx)*
                 up.mobility(phaseIdx)*
                 up.enthalpy(phaseIdx)
                 +
                 (1 - mobilityUpwindAlpha)*
                 dn.density(phaseIdx)*
                 dn.mobility(phaseIdx)*
                 dn.enthalpy(phaseIdx));
        }
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxData The flux variables at the current SCV
     *
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxData) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxData);

        // diffusive heat flux
        flux[energyEqIdx] += fluxData.normalMatrixHeatFlux();
    }
};

}

#endif
