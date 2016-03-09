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
 *        using the one-phase fully implicit model.
 */
#ifndef DUMUX_1P_LOCAL_RESIDUAL_HH
#define DUMUX_1P_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup OnePModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase fully implicit model.
 */
template<class TypeTag>
class OnePLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    //index of the mass balance equation
    enum {
        conti0EqIdx = Indices::conti0EqIdx //index for the mass balance
    };

public:

    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    OnePLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the OneP
     *        model.
     *
     * This function should not include the source and sink terms.
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, const bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // partial time derivative of the wetting phase mass
        storage[conti0EqIdx] =  volVars.density() * volVars.porosity();
    }


    /*!
     * \brief Evaluate the mass flux over a face of a sub-control
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param fIdx The index of the SCV face
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

        asImp_()->computeAdvectiveFlux(flux, fluxVars);
        asImp_()->computeDiffusiveFlux(flux, fluxVars);
    }

    /*!
     * \brief Evaluate the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(/*phaseIdx=*/0));
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(/*phaseIdx=*/0));
        flux[conti0EqIdx] =
            ((    upwindWeight_)*up.density()
             +
             (1 - upwindWeight_)*dn.density())
            *
            fluxVars.volumeFlux(/*phaseIdx=*/0);
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // diffusive fluxes
        flux += 0.0;
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class PrimaryVariables>
    Scalar temperature(const PrimaryVariables &priVars)
    { return this->problem_.temperature(); /* constant temperature */ }

private:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

    Scalar upwindWeight_;
};

}

#endif
