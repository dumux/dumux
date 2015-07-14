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
 * \brief Element-wise calculation of the residual for the Richards fully implicit model.
 */
#ifndef DUMUX_RICHARDS_LOCAL_RESIDUAL_HH
#define DUMUX_RICHARDS_LOCAL_RESIDUAL_HH

#include "richardsproperties.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the residual for the Richards fully implicit model.
 */
template<class TypeTag>
class RichardsLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        contiEqIdx = Indices::contiEqIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    static const bool useHead = GET_PROP_VALUE(TypeTag, UseHead);

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    RichardsLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the Richards
     *        model.
     *
     * This function should not include the source and sink terms.
     *
     * \param storage Stores the average mass per unit volume for each phase \f$\mathrm{[kg/m^3]}\f$
     * \param scvIdx The sub control volume index of the current element
     * \param usePrevSol Calculate the storage term of the previous solution
     *                   instead of the model's current solution.
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VolumeVariables &volVars =
            usePrevSol ?
            this->prevVolVars_(scvIdx) :
            this->curVolVars_(scvIdx);

        // partial time derivative of the wetting phase mass
        //pressure head formulation
        storage[contiEqIdx] =
                volVars.saturation(wPhaseIdx)
                * volVars.porosity();
        // pressure formulation
        if (!useHead)
            storage[contiEqIdx] *= volVars.density(wPhaseIdx);
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     *
     * \param flux Stores the total mass fluxes over a sub-control volume face
     *             of the current element \f$\mathrm{[kg/s]}\f$
     * \param fIdx The sub control volume face index inside the current
     *                element
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, bool onBoundary=false) const
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
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */

    void computeAdvectiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(wPhaseIdx));
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(wPhaseIdx));

        //pressure head formulation
        flux[contiEqIdx] =
            fluxVars.volumeFlux(wPhaseIdx);
        //pressure formulation
        if(!useHead)
            flux[contiEqIdx] *=((    massUpwindWeight_)*up.density(wPhaseIdx)
                                +
                                (1 - massUpwindWeight_)*dn.density(wPhaseIdx));
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * This function doesn't do anything but may be used by the
     * non-isothermal three-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(PrimaryVariables &flux, const FluxVariables &fluxVars) const
    {
        // diffusive fluxes
        flux += 0.0;
    }


    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source Stores the average source term of all phases inside a
     *          sub-control volume of the current element \f$\mathrm{[kg/(m^3 * s)]}\f$
     * \param scvIdx The sub control volume index inside the current
     *               element
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

private:
    Scalar massUpwindWeight_;


};

}

#endif
