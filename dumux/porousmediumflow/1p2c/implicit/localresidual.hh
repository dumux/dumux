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
 * \brief Element-wise calculation the local Jacobian for the single-phase,
 *        two-component model in the fully implicit scheme.
 */

#ifndef DUMUX_ONEP_TWOC_LOCAL_RESIDUAL_HH
#define DUMUX_ONEP_TWOC_LOCAL_RESIDUAL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 *
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitLocalResidual
 * \brief Calculate the local Jacobian for the single-phase,
 *        two-component model in the fully implicit scheme.
 *
 *  This class is used to fill the gaps in BoxLocalResidual for the 1p2c flow and transport.
 */
template<class TypeTag>
class OnePTwoCLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        //phase index
        phaseIdx = Indices::phaseIdx,
        transportCompIdx = Indices::transportCompIdx
    };
    // indices of the equations
    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx
    };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);



public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    OnePTwoCLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a finite volume.
     *
     *        \param storage The mass of the component within the sub-control volume
     *        \param scvIdx The index of the considered face of the sub-control volume
     *        \param usePrevSol Evaluate function with solution of current or previous time step
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

        storage = 0;
        if(useMoles) // mole-fraction formulation
        {
            // storage term of continuity equation- molefractions
            //careful: molarDensity changes with moleFrac!
            storage[conti0EqIdx] += volVars.molarDensity()*volVars.porosity();
            // storage term of the transport equation - molefractions
            storage[transportEqIdx] +=
                volVars.molarDensity()*volVars.moleFraction(transportCompIdx) *
                volVars.porosity();
        }
        else // mass-fraction formulation
        {
            // storage term of continuity equation - massfractions
            storage[conti0EqIdx] +=
                volVars.density()*volVars.porosity();
            //storage term of the transport equation - massfractions
            storage[transportEqIdx] +=
                volVars.density() * volVars.massFraction(transportCompIdx) * volVars.porosity();
        }
    }

    /*!
     * \brief Evaluate the mass flux over a face of a sub-control
     *        volume.
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param fIdx The index of the considered face of the sub control volume
     *        \param onBoundary A boolean variable to specify whether the flux variables
     *               are calculated for interior SCV faces or boundary faces, default=false
     */
    void computeFlux(PrimaryVariables &flux, const int fIdx, const bool onBoundary=false) const
    {
        flux = 0;
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
        ////////
        // advective fluxes of all components in all phases
        ////////

        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up =
            this->curVolVars_(fluxVars.upstreamIdx());
        const VolumeVariables &dn =
            this->curVolVars_(fluxVars.downstreamIdx());

        if(!useMoles) //mass-fraction formulation
        {
            // total mass flux - massfraction
            //KmvpNormal is the Darcy velocity multiplied with the normal vector, calculated in 1p2cfluxvariables.hh
            flux[conti0EqIdx] +=
                fluxVars.KmvpNormal() *
                ((     upwindWeight_)*up.density()/up.viscosity()
                 +
                 ((1 - upwindWeight_)*dn.density()/dn.viscosity()));

            // advective flux of the second component - massfraction
            flux[transportEqIdx] +=
                fluxVars.KmvpNormal() *
                ((    upwindWeight_)*up.density() * up.massFraction(transportCompIdx)/up.viscosity()
                 +
                 (1 - upwindWeight_)*dn.density()*dn.massFraction(transportCompIdx)/dn.viscosity());
        }
        else //mole-fraction formulation
        {
            // total mass flux - molefraction
            //KmvpNormal is the Darcy velocity multiplied with the normal vector, calculated in 1p2cfluxvariables.hh
            flux[conti0EqIdx] +=
                fluxVars.KmvpNormal() *
                ((     upwindWeight_)*up.molarDensity()/up.viscosity()
                 +
                 ((1 - upwindWeight_)*dn.molarDensity()/dn.viscosity()));

            // advective flux of the second component -molefraction
            flux[transportEqIdx] +=
                fluxVars.KmvpNormal() *
                ((    upwindWeight_)*up.molarDensity() * up.moleFraction(transportCompIdx)/up.viscosity()
                 +
                 (1 - upwindWeight_)*dn.molarDensity() * dn.moleFraction(transportCompIdx)/dn.viscosity());
        }

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
        Scalar tmp(0);

        // diffusive flux of second component
        if(useMoles) // mole-fraction formulation
        {
            // diffusive flux of the second component - molefraction
            tmp = -(fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal);
            tmp *= fluxVars.porousDiffCoeff() * fluxVars.molarDensity();

            // dispersive flux of second component - molefraction
                        GlobalPosition normalDisp;
                        fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
                        tmp -= fluxVars.molarDensity()*
                            (normalDisp * fluxVars.moleFractionGrad(transportCompIdx));

            flux[transportEqIdx] += tmp;
        }
        else // mass-fraction formulation
        {
            // diffusive flux of the second component - massfraction
            tmp = -(fluxVars.moleFractionGrad(transportCompIdx)*fluxVars.face().normal);
            tmp *= fluxVars.porousDiffCoeff() * fluxVars.molarDensity();

            // dispersive flux of second component - massfraction
                       GlobalPosition normalDisp;
                       fluxVars.dispersionTensor().mv(fluxVars.face().normal, normalDisp);
                       tmp -= fluxVars.molarDensity()*
                       (normalDisp * fluxVars.moleFractionGrad(transportCompIdx));

            // convert it to a mass flux and add it
            flux[transportEqIdx] += tmp * FluidSystem::molarMass(transportCompIdx);
        }
    }

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }

private:
    Scalar upwindWeight_;
};

}

#endif
