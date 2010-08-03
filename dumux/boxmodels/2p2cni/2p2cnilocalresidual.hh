// $Id: 2p2cnilocalresidual.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal two-phase two-component box model.
 *
 */
#ifndef DUMUX_NEW_2P2CNI_BOX_JACOBIAN_HH
#define DUMUX_NEW_2P2CNI_BOX_JACOBIAN_HH

#include <dumux/boxmodels/2p2c/2p2clocalresidual.hh>


#include <dumux/boxmodels/2p2cni/2p2cnisecondaryvars.hh>
#include <dumux/boxmodels/2p2cni/2p2cnifluxvars.hh>

#include <dumux/boxmodels/2p2cni/2p2cniproperties.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCNIModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component box model.
 */
template<class TypeTag>
class TwoPTwoCNILocalResidual : public TwoPTwoCLocalResidual<TypeTag>
{
    typedef TwoPTwoCNILocalResidual<TypeTag> ThisType;
    typedef TwoPTwoCLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVarVector)) PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    enum {
        dim              = GridView::dimension,
        dimWorld         = GridView::dimensionworld,

        numPhases        = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        temperatureIdx   = Indices::temperatureIdx,

        lPhaseIdx   = Indices::lPhaseIdx,
        gPhaseIdx   = Indices::gPhaseIdx
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) SecondaryVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVars)) FluxVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSecondaryVars)) ElementSecondaryVars;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    static const Scalar mobilityUpwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(PrimaryVarVector &result, int scvIdx, bool usePrevSol) const
    {
        // compute the storage term for phase mass
        ParentType::computeStorage(result, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementSecondaryVars &elemDat = usePrevSol ? this->prevSecVars_()  : this->curSecVars_();
        const SecondaryVars  &vertDat = elemDat[scvIdx];

        // compute the energy storage
        result[temperatureIdx] =
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
            vertDat.temperature() *
            vertDat.heatCapacity();
    }

    /*!
     * \brief Evaluates the advective mass flux and the heat flux
     * over a face of a subcontrol volume and writes the result in
     * the flux vector.
     *
     * This method is called by compute flux (base class)
     */
    void computeAdvectiveFlux(PrimaryVarVector &flux,
                              const FluxVars &fluxData) const
    {
        // advective mass flux
        ParentType::computeAdvectiveFlux(flux, fluxData);

        // advective heat flux in all phases
        flux[temperatureIdx] = 0;
        for (int phase = 0; phase < numPhases; ++phase) {
            // vertex data of the upstream and the downstream vertices
            const SecondaryVars &up = this->curSecVars_(fluxData.upstreamIdx(phase));
            const SecondaryVars &dn = this->curSecVars_(fluxData.downstreamIdx(phase));

            flux[temperatureIdx] +=
                fluxData.KmvpNormal(phase) * (
                    mobilityUpwindAlpha * // upstream vertex
                    (  up.density(phase) *
                       up.mobility(phase) *
                       up.enthalpy(phase))
                    +
                    (1-mobilityUpwindAlpha) * // downstream vertex
                    (  dn.density(phase) *
                       dn.mobility(phase) *
                       dn.enthalpy(phase)) );
        }
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     */
    void computeDiffusiveFlux(PrimaryVarVector &flux,
                              const FluxVars &fluxData) const
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
