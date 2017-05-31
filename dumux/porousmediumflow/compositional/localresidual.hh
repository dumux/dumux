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
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */
#ifndef DUMUX_COMPOSITIONAL_LOCAL_RESIDUAL_HH
#define DUMUX_COMPOSITIONAL_LOCAL_RESIDUAL_HH

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(ReplaceCompEqIdx);
} // end namespace Properties

/*!
 * \ingroup Implicit
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 *
 */
template<class TypeTag>
class CompositionalLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseLocalResidual);
    using Implementation = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using EnergyLocalResidual = typename GET_PROP_TYPE(TypeTag, EnergyLocalResidual);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    enum { conti0EqIdx = Indices::conti0EqIdx };

    //! The index of the component balance equation that gets replaced with the total mass balance
    static const int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);

public:

    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    PrimaryVariables computeStorage(const SubControlVolume& scv,
                                    const VolumeVariables& volVars) const
    {
        PrimaryVariables storage(0.0);

        // formulation with mole balances
        if (useMoles)
        {
            // compute storage term of all components within all phases
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    auto eqIdx = conti0EqIdx + compIdx;
                    if (eqIdx != replaceCompEqIdx)
                        storage[eqIdx] += volVars.porosity()
                                          * volVars.saturation(phaseIdx)
                                          * volVars.molarDensity(phaseIdx)
                                          * volVars.moleFraction(phaseIdx, compIdx);
                }

                // in case one balance is substituted by the total mole balance
                if (replaceCompEqIdx < numComponents)
                    storage[replaceCompEqIdx] = volVars.molarDensity(phaseIdx)
                                                * volVars.porosity()
                                                * volVars.saturation(phaseIdx);

                //! The energy storage in the fluid phase with index phaseIdx
                EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
            }

            //! The energy storage in the solid matrix
            EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);
        }
        // formulation with mass balances
        else
        {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    auto eqIdx = conti0EqIdx + compIdx;
                    if (eqIdx != replaceCompEqIdx)
                        storage[eqIdx] += volVars.porosity()
                                          * volVars.saturation(phaseIdx)
                                          * volVars.density(phaseIdx)
                                          * volVars.massFraction(phaseIdx, compIdx);
                }

                // in case one balance is substituted by the total mass balance
                if (replaceCompEqIdx < numComponents)
                    storage[replaceCompEqIdx] = volVars.density(phaseIdx)
                                                * volVars.porosity()
                                                * volVars.saturation(phaseIdx);

                //! The energy storage in the fluid phase with index phaseIdx
                EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
            }

            //! The energy storage in the solid matrix
            EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);
        }

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each component
     * \param fIdx The index of the SCV face
     * \param onBoundary A boolean variable to specify whether the flux variables
     *        are calculated for interior SCV faces or boundary faces, default=false
     */
    PrimaryVariables computeFlux(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        FluxVariables fluxVars;
        fluxVars.init(this->problem(), element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        PrimaryVariables flux(0.0);

        // formulation with mole balances
        if (useMoles)
        {
            // advective fluxes
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    // get equation index
                    auto eqIdx = conti0EqIdx + compIdx;

                    // the physical quantities for which we perform upwinding
                    auto upwindTerm = [phaseIdx, compIdx](const auto& volVars)
                    { return volVars.molarDensity(phaseIdx)*volVars.moleFraction(phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                    if (eqIdx != replaceCompEqIdx)
                        flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                    // diffusive fluxes (only for the component balances)
                    if (phaseIdx != compIdx)
                    {
                        //one of these if is always true
                        if (eqIdx != replaceCompEqIdx)
                            flux[eqIdx] += diffusiveFluxes[eqIdx];
                        if (conti0EqIdx + phaseIdx != replaceCompEqIdx)
                            flux[conti0EqIdx + phaseIdx] -= diffusiveFluxes[conti0EqIdx+ phaseIdx];
                    }
                }

                // in case one balance is substituted by the total mole balance
                if (replaceCompEqIdx < numComponents)
                {
                    auto upwindTermTotalBalance = [phaseIdx](const auto& volVars)
                    { return volVars.molarDensity(phaseIdx)*volVars.mobility(phaseIdx); };

                    flux[replaceCompEqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindTermTotalBalance);
                }

                //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
                EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
            }

            //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConductionFlux(flux, fluxVars);
        }
        // formulation with mass balances
        else
        {
            // advective fluxes
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                const auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    // get equation index
                    auto eqIdx = conti0EqIdx + compIdx;

                    // the physical quantities for which we perform upwinding
                    auto upwindTerm = [phaseIdx, compIdx](const auto& volVars)
                    { return volVars.density(phaseIdx)*volVars.massFraction(phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                    const auto advFlux = fluxVars.advectiveFlux(phaseIdx, upwindTerm);

                    if (eqIdx != replaceCompEqIdx)
                        flux[eqIdx] += advFlux;

                    // diffusive fluxes (only for the component balances)
                    if (phaseIdx != compIdx)
                    {
                        if (eqIdx != replaceCompEqIdx)
                            flux[eqIdx] +=  diffusiveFluxes[eqIdx]*FluidSystem::molarMass(compIdx);
                        if (conti0EqIdx + phaseIdx != replaceCompEqIdx)
                            flux[conti0EqIdx + phaseIdx] -=  diffusiveFluxes[conti0EqIdx+ phaseIdx]*FluidSystem::molarMass(compIdx);
                    }
                }

                // in case one balance is substituted by the total mass balance
                if (replaceCompEqIdx < numComponents)
                {
                    auto upwindTermTotalBalance = [phaseIdx](const auto& volVars)
                    { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx); };

                    flux[replaceCompEqIdx] = fluxVars.advectiveFlux(phaseIdx, upwindTermTotalBalance);

                }

                //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
                EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
            }

            //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConductionFlux(flux, fluxVars);
        }

        return flux;
    }

protected:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

} // end namespace Dumux

#endif
