// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 *        using the compositional stokes box model.
 *
 */
#ifndef DUMUX_STOKES2C_LOCAL_RESIDUAL_HH
#define DUMUX_STOKES2C_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/stokes/stokeslocalresidual.hh>

#include <dumux/freeflow/stokes2c/stokes2cvolumevariables.hh>
#include <dumux/freeflow/stokes2c/stokes2cfluxvariables.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the compositional stokes box model. This is derived
 *        from the stokes box model.
 */
template<class TypeTag>
class Stokes2cLocalResidual : public StokesLocalResidual<TypeTag>
{
    typedef StokesLocalResidual<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cIndices) Indices;

    enum { dim = GridView::dimension };
    enum {
        transportIdx = Indices::transportIdx //!< Index of the transport equation
    };
    enum {
        lCompIdx = Indices::lCompIdx
    };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIndex)};

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;


public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    Stokes2cLocalResidual()
    {};

    /*!
     * \brief Evaluate the amount the additional quantities to the stokes model
     *        (transport equation).
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(PrimaryVariables &result, int scvIdx, bool usePrevSol) const
    {
        // compute the storage term for the transport equation
        ParentType::computeStorage(result, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemDat = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &vertexDat = elemDat[scvIdx];

        // compute the storage of the component
        result[transportIdx] =
            vertexDat.density() *
            vertexDat.fluidState().massFraction(phaseIdx, lCompIdx);

        Valgrind::CheckDefined(vertexDat.density());
        Valgrind::CheckDefined(vertexDat.fluidState().massFraction(phaseIdx, lCompIdx));
    }

    /*!
     * \brief Evaluates the advective component (mass) flux
     * over a face of a subcontrol volume and writes the result in
     * the flux vector.
     *
     * This method is called by compute flux (base class)
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // call computation of the advective fluxes of the stokes model
        // (momentum and mass fluxes)
        ParentType::computeAdvectiveFlux(flux, fluxVars);

        // vertex data of the upstream and the downstream vertices
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx());
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx());

        Scalar tmp = fluxVars.normalVelocityAtIP();

        if (this->massUpwindWeight_ > 0.0)
            tmp *=  this->massUpwindWeight_ *         // upwind data
                up.density() * up.fluidState().massFraction(phaseIdx, lCompIdx);
        if (this->massUpwindWeight_ < 1.0)
            tmp += (1.0 - this->massUpwindWeight_) *     // rest
                dn.density() * dn.fluidState().massFraction(phaseIdx, lCompIdx);

        flux[transportIdx] += tmp;
        Valgrind::CheckDefined(flux[transportIdx]);
    }

    /*!
     * \brief Adds the diffusive component flux to the flux vector over
     *        the face of a sub-control volume.
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxVars);

        // diffusive component flux
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            flux[transportIdx] -=
                fluxVars.massFractionGradAtIP()[dimIdx] *
                fluxVars.face().normal[dimIdx] *
                fluxVars.diffusionCoeffAtIP() *
                fluxVars.densityAtIP();
        Valgrind::CheckDefined(flux[transportIdx]);
    }
};

}

#endif
