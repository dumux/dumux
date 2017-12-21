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
 * \ingroup RichardsModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the Richards fully implicit models.
 */
#ifndef DUMUX_RICHARDS_LOCAL_RESIDUAL_HH
#define DUMUX_RICHARDS_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the Richards fully implicit models.
 */
template<class TypeTag>
class RichardsLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);

    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
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
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    // first index for the mass balance
    enum { conti0EqIdx = Indices::conti0EqIdx };

    // phase indices
    enum {
           wPhaseIdx = FluidSystem::wPhaseIdx,
           nPhaseIdx = FluidSystem::nPhaseIdx,
           wCompIdx = FluidSystem::wCompIdx
    };

    static constexpr bool enableWaterDiffusionInAir
        = GET_PROP_VALUE(TypeTag, EnableWaterDiffusionInAir);
public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    ResidualVector computeStorage(const Problem& problem,
                                  const SubControlVolume& scv,
                                  const VolumeVariables& volVars) const
    {
        // partial time derivative of the phase mass
        ResidualVector storage(0.0);
        storage[conti0EqIdx] = volVars.porosity()
                               * volVars.density(wPhaseIdx)
                               * volVars.saturation(wPhaseIdx);

        // for extended Richards we consider water in air
        if (enableWaterDiffusionInAir)
            storage[conti0EqIdx] += volVars.porosity()
                                    * volVars.molarDensity(nPhaseIdx)
                                    * volVars.moleFraction(nPhaseIdx, wCompIdx)
                                    * FluidSystem::molarMass(wCompIdx)
                                    * volVars.saturation(nPhaseIdx);

        //! The energy storage in the water, air and solid phase
        EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, wPhaseIdx);
        EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, nPhaseIdx);
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }


    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux compuation
     */
    ResidualVector computeFlux(const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        ResidualVector flux(0.0);
        // the physical quantities for which we perform upwinding
        auto upwindTerm = [](const auto& volVars)
                          { return volVars.density(wPhaseIdx)*volVars.mobility(wPhaseIdx); };

        flux[conti0EqIdx] = fluxVars.advectiveFlux(wPhaseIdx, upwindTerm);

        // for extended Richards we consider water vapor diffusion in air
        if (enableWaterDiffusionInAir)
            flux[conti0EqIdx] += fluxVars.molecularDiffusionFlux(nPhaseIdx)[wCompIdx]*FluidSystem::molarMass(wCompIdx);

        //! Add advective phase energy fluxes for the water phase only. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, wPhaseIdx);

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        //! The effective lambda is averaged over both fluid phases and the solid phase
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
