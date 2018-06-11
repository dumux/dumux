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
 *        using the n-phase immiscible fully implicit models.
 */
#ifndef DUMUX_TWOP_FRACFLOW_LOCAL_RESIDUAL_HH
#define DUMUX_TWOP_FRACFLOW_LOCAL_RESIDUAL_HH

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the n-phase immiscible fully implicit models.
 */
template<class TypeTag>
class TwoPFractionalFlowLocalResidual
: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    // using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);

    // first index for the mass balance
    enum {
           transportEqIdx = 0,
           totalvelocityEqIdx = 1,
           wPhaseIdx = 0,
           nPhaseIdx = 1
    };

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    template<class Problem, class SubControlVolume, class VolumeVariables>
    ResidualVector computeStorage(const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars) const
    {
        // partial time derivative of the phase mass
        ResidualVector storage;

        // the saturation transport equation storage term
        storage[transportEqIdx] = volVars.porosity()
                                  // * volVars.density(wPhaseIdx)
                                  * volVars.saturation(wPhaseIdx);
        storage[totalvelocityEqIdx] = 0;

        //! The energy storage in the fluid phase with index phaseIdx
        // EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, wPhaseIdx);
        // EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, nPhaseIdx);

        //! The energy storage in the solid matrix
        // EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }


    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume
     * \param scvf The sub control volume face to compute the flux on
     */
    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class ElementFluxVariablesCache>
    ResidualVector computeFlux(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        ResidualVector flux;

        // this is just a dummy we do upwinding internally in the upwind scheme class
        auto upwindTermTransport = [](const auto& volVars, const int phaseIdx = 0){ return 0.0; };
        flux[transportEqIdx] = fluxVars.advectiveFlux(transportEqIdx, upwindTermTransport);
        flux[totalvelocityEqIdx] = fluxVars.advetiveFlux(totalvelocityEqIdx,upwindTermTransport) +fluxVars.advectiveFlux(transportEqIdx, upwindTermTransport);

        //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
        // EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, wPhaseIdx);
        // EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, nPhaseIdx);

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        // EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;
    }
};

} // end namespace Dumux

#endif
