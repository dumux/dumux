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
 *        using the non-isothermal Stokes box model with a staggered grid.
 *
 */
#ifndef DUMUX_STAGGERED_NAVIERSTOKES_NI_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_NAVIERSTOKES_NI_LOCAL_RESIDUAL_HH

//#include <dune/istl/matrix.hh> // TODO ?

//#include <dumux/common/valgrind.hh> // TODO ?
#include <dumux/implicit/staggered/localresidual.hh>

//#include "properties.hh" // TODO ?
//
//#include "volumevariables.hh" // TODO ?
//#include "fluxvariables.hh" // TODO ?

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableEnergyBalance);
NEW_PROP_TAG(EnableInertiaTerms);
//NEW_PROP_TAG(ReplaceCompEqIdx);
}

/*!
 * \ingroup CCModel
 * \ingroup StaggeredLocalResidual
 * \brief Element-wise calculation of the residual for models
 * 			based on the fully implicit cell-centered scheme
 *
 * \todo Please doc me more!
 */

// forward declaration
template<class TypeTag, bool enableComponentTransport, bool enableEnergyBalance>
class StaggeredNavierStokesResidualImpl;

// specialization for immiscible, non-isothermal flow
template<class TypeTag>
class StaggeredNavierStokesResidualImpl<TypeTag, false, true> : public StaggeredNavierStokesResidualImpl<TypeTag, false, false>
{
    friend class StaggeredLocalResidual<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum { // TODO adapt
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,
        velocityIdx = Indices::velocityIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        energyBalanceIdx = Indices::energyBalanceIdx
    };

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);

//    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief Evaluate the amount the additional quantities to the stokes model
     *        (energy equation).
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    CellCenterPrimaryVariables computeStorageForCellCenter(const SubControlVolume& scv, const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage(0.0);

//        const Scalar density = useMoles? volVars.molarDensity() : volVars.density();

        // compute storage of mass
        storage[massBalanceIdx] = volVars.density(0);

        // compute the storage of energy
        storage[energyBalanceIdx] = volVars.density(0) * volVars.internalEnergy(0);

        return storage;
    }

    // TODO implement advectiveFlux, conductiveFlux

    /*!
     * \brief Evaluates the convective energy flux
     *        over a face of a sub-control volume and writes the result in
     *        the flux vector. This method is called by computeFlux in the base class.
     *
     * \param flux The vector for the fluxes over the SCV/boundary face
     * \param fluxVars The flux variables at the current SCV/boundary face
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

        Scalar tmp = fluxVars.normalVelocity();
        tmp *= (this->massUpwindWeight_ * up.density() * up.enthalpy()
                + (1.0 - this->massUpwindWeight_) * dn.density() * dn.enthalpy());

        flux[energyEqIdx] += tmp;
        Valgrind::CheckDefined(flux[energyEqIdx]);
    }

    /*!
     * \brief Evaluates the diffusive component energy flux
     *        over the face of a sub-control volume.
     *
     * \param flux The vector for the fluxes over the SCV face
     * \param fluxVars The flux variables at the current SCV face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxVars);

        // conductive energy flux
        computeConductiveFlux(flux, fluxVars);

        // diffusive component energy flux
        Scalar sumDiffusiveFluxes = 0;
        for (int compIdx=0; compIdx<numComponents; compIdx++)
        {
            if (compIdx != phaseCompIdx)
            {
                Valgrind::CheckDefined(fluxVars.moleFractionGrad(compIdx));
                Valgrind::CheckDefined(fluxVars.face().normal);
                Valgrind::CheckDefined(fluxVars.diffusionCoeff(compIdx));
                Valgrind::CheckDefined(fluxVars.eddyDiffusivity());
                Valgrind::CheckDefined(fluxVars.molarDensity());
                Valgrind::CheckDefined(FluidSystem::molarMass(compIdx));
                Valgrind::CheckDefined(fluxVars.componentEnthalpy(compIdx));
                Scalar diffusiveFlux = fluxVars.moleFractionGrad(compIdx)
                  * fluxVars.face().normal
                  * (fluxVars.diffusionCoeff(compIdx) + fluxVars.eddyDiffusivity())
                  * fluxVars.molarDensity();
                sumDiffusiveFluxes += diffusiveFlux;
                flux[energyEqIdx] -= diffusiveFlux * fluxVars.componentEnthalpy(compIdx)
                  * FluidSystem::molarMass(compIdx); // Multiplied by molarMass [kg/mol] to convert from [mol/m^3 s] to [kg/m^3 s];
            }
        }

        // the diffusive flux of the phase component is the negative of the sum of the component fluxes
        flux[energyEqIdx] += sumDiffusiveFluxes * fluxVars.componentEnthalpy(phaseCompIdx)
          * FluidSystem::molarMass(phaseCompIdx); // Multiplied by molarMass [kg/mol] to convert from [mol/m^3 s] to [kg/m^3 s];

        Valgrind::CheckDefined(flux[energyEqIdx]);
    }

    /*!
     * \brief Evaluates the conductive energy flux
     *        over the face of a sub-control volume.
     *
     * \param flux The vector for the fluxes over the SCV face
     * \param fluxVars The flux variables at the current SCV face
     */
    void computeConductiveFlux(PrimaryVariables &flux,
                               const FluxVariables &fluxVars) const
    {
        // diffusive heat flux
        flux[energyEqIdx] -=
            fluxVars.temperatureGrad() * fluxVars.face().normal
            * (fluxVars.thermalConductivity() + fluxVars.thermalEddyConductivity());
        Valgrind::CheckDefined(flux[energyEqIdx]);
    }
};

} // end namespace

#endif
