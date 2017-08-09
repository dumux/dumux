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

namespace Dumux
{
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the n-phase immiscible fully implicit models.
 */
template<class TypeTag>
class TwoPFractionalFlowLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // first index for the mass balance
    enum { conti0EqIdx = Indices::conti0EqIdx,
           wPhaseIdx = Indices::wPhaseIdx };

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

public:

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
    PrimaryVariables computeStorage(const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        // partial time derivative of the phase mass
        PrimaryVariables storage;

        storage[conti0EqIdx] = volVars.porosity()
                                * volVars.density(wPhaseIdx)
                                * volVars.saturation(wPhaseIdx);

        //! The energy storage in the fluid phase with index phaseIdx
        EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, wPhaseIdx);

        //! The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }


    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume
     * \param scvf The sub control volume face to compute the flux on
     */
    PrimaryVariables computeFlux(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(this->problem(), element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        PrimaryVariables flux;

        // this is just a dummy we do upwinding internally in the upwind scheme class
        auto upwindTerm = [](const auto& volVars){};

        flux[conti0EqIdx] = fluxVars.advectiveFlux(wPhaseIdx, upwindTerm);

        //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, wPhaseIdx);

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;
    }

private:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

} // end namespace Dumux

#endif
