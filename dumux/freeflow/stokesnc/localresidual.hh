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
 *        using the compositional Stokes box model.
 *
 */
#ifndef DUMUX_STOKESNC_LOCAL_RESIDUAL_HH
#define DUMUX_STOKESNC_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/stokes/localresidual.hh>

#include "volumevariables.hh"
#include "fluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup BoxStokesncModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the compositional Stokes box model. This is derived
 *        from the Stokes box model.
 */
template<class TypeTag>
class StokesncLocalResidual : public StokesLocalResidual<TypeTag>
{
    typedef StokesLocalResidual<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    //dimensions
    enum {  dim = GridView::dimension };
    //number of equations
    enum {  numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    //number of components
    enum {  numComponents = Indices::numComponents };
    //equation indices
    enum {  massBalanceIdx = Indices::massBalanceIdx,
            momentumXIdx = Indices::momentumXIdx,
            lastMomentumIdx = Indices::lastMomentumIdx,
            transportEqIdx = Indices::transportEqIdx,
            conti0EqIdx = Indices::conti0EqIdx };
    //primary variable indices
    enum {  pressureIdx = Indices::pressureIdx };
    //phase employed
    enum {  phaseIdx = Indices::phaseIdx };
    //component indices
    enum {  phaseCompIdx = Indices::phaseCompIdx,
            transportCompIdx = Indices::transportCompIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    static const bool calculateNavierStokes = GET_PROP_VALUE(TypeTag, EnableNavierStokes);

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief Evaluate the stored amount of quantities additional to the Stokes model
     *        (transport equations). For using mole fraction also momentum balances and mass balance
     *         have to be calculated using molar quantities
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, const bool usePrevSol) const
    {
        // compute the storage term for the transport equation
        ParentType::computeStorage(storage, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit Euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ?
                    this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        if (useMoles)
        {
            // mass balance and transport equations
            for (int compIdx=0; compIdx<numComponents; compIdx++)
            {
                if (conti0EqIdx+compIdx != massBalanceIdx)
                //else // transport equations
                {
                    storage[conti0EqIdx+compIdx] = volVars.molarDensity()
                        * volVars.moleFraction(compIdx);

                    Valgrind::CheckDefined(volVars.molarDensity());
                    Valgrind::CheckDefined(volVars.moleFraction(compIdx));
                }
            }
        }
        else
        {
            /* works for a maximum of two components, for more components
            mole fractions must be used (set property UseMoles to true) */

            //storage of transported component
            storage[transportEqIdx] = volVars.density()
                * volVars.massFraction(transportCompIdx);

            Valgrind::CheckDefined(volVars.density());
            Valgrind::CheckDefined(volVars.massFraction(transportCompIdx));

        }
    }

    /*!
     * \brief Evaluates the advective component (mass) flux
     * over a face of a sub-control volume and writes the result in
     * the flux vector.
     *
     * This method is called by compute flux (base class).
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // call ParentType function
        ParentType::computeAdvectiveFlux(flux,fluxVars);

        // data attached to upstream and the downstream vertices
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx());
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx());

        Scalar tmp = fluxVars.normalVelocity();

        if(useMoles)
        {
            //transport equations
            for (int compIdx=0; compIdx<numComponents; compIdx++)
            {
                if (conti0EqIdx+compIdx != massBalanceIdx) // mass balance is calculated above
                {
                    tmp *= (this->massUpwindWeight_ * up.molarDensity() * up.moleFraction(compIdx)
                            + (1.-this->massUpwindWeight_) * dn.molarDensity() * dn.moleFraction(compIdx));

                    flux[conti0EqIdx+compIdx] += tmp;
                    Valgrind::CheckDefined(flux[conti0EqIdx+compIdx]);
                }
            }
        }
        else
        {
            tmp *= (this->massUpwindWeight_ * up.density() * up.massFraction(transportCompIdx)
                    + (1.-this->massUpwindWeight_) * dn.density() * dn.massFraction(transportCompIdx));

            flux[transportEqIdx] += tmp;
            Valgrind::CheckDefined(flux[transportEqIdx]);
        }
    }

    /*!
     * \brief Adds the diffusive component flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the SCV face or boundary face
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // diffusive component flux
        if(useMoles)
        {
            //loop over secondary components
            for (int compIdx=0; compIdx<numComponents; compIdx++)
            {
                if (conti0EqIdx+compIdx != massBalanceIdx)
                {
                    flux[conti0EqIdx+compIdx] += -(fluxVars.moleFractionGrad(compIdx) * fluxVars.face().normal)
                                                 * (fluxVars.diffusionCoeff(compIdx) + fluxVars.eddyDiffusivity())
                                                 * fluxVars.molarDensity();
                    Valgrind::CheckDefined(flux[conti0EqIdx+compIdx]);
                }
            }
        }
        else
        {
            Scalar diffusiveMolarFlux = -(fluxVars.moleFractionGrad(transportCompIdx) * fluxVars.face().normal)
                                        * (fluxVars.diffusionCoeff(transportCompIdx) + fluxVars.eddyDiffusivity())
                                        * fluxVars.molarDensity();
            // Multiply by molarMass [kg/mol] to convert form [mol/m^3 s] to [kg/m^3 s]
            flux[transportEqIdx] += diffusiveMolarFlux * FluidSystem::molarMass(transportCompIdx);
            // Add the diffusive fluxes to the total mass balance, this is necessary as
            // diffusive mole fluxes cancel out, but not the diffusive mass fluxes
            // NOTE: for the phaseCompIdx the diffusiveMolarFlux is inverted, because it was calculated
            //       by moleFractionGrad(transportCompIdx)
            flux[massBalanceIdx] += diffusiveMolarFlux * FluidSystem::molarMass(transportCompIdx);
            flux[massBalanceIdx] += -diffusiveMolarFlux * FluidSystem::molarMass(phaseCompIdx);
            Valgrind::CheckDefined(flux[transportEqIdx]);
        }
    }
};

}

#endif
