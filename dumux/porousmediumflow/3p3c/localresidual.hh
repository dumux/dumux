// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup ThreePThreeCModel
  * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase three-component fully implicit model.
 */

#ifndef DUMUX_3P3C_LOCAL_RESIDUAL_HH
#define DUMUX_3P3C_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase three-component fully implicit model.
 */
template<class TypeTag>
class ThreePThreeCLocalResidual : public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum {
        numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases(),
        numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents(),

        contiWEqIdx = Indices::conti0EqIdx + FluidSystem::wPhaseIdx,//!< index of the mass conservation equation for the water component
        contiNEqIdx = Indices::conti0EqIdx + FluidSystem::nPhaseIdx,//!< index of the mass conservation equation for the contaminant component
        contiGEqIdx = Indices::conti0EqIdx + FluidSystem::gPhaseIdx,//!< index of the mass conservation equation for the gas component

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,
        gCompIdx = FluidSystem::gCompIdx
    };

public:

    using ParentType::ParentType;

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The volume variables
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                auto eqIdx = Indices::conti0EqIdx + compIdx;
                storage[eqIdx] += volVars.porosity()
                                  * volVars.saturation(phaseIdx)
                                  * volVars.molarDensity(phaseIdx)
                                  * volVars.moleFraction(phaseIdx, compIdx);
            }

            // The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        // The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume element geometry
     * \param elemVolVars The element volume variables
     * \param scvf The sub control volume face
     * \param elemFluxVarsCache The element flux variables cache
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        static constexpr auto referenceSystemFormulation = FluxVariables::MolecularDiffusionType::referenceSystemFormulation();

        // get upwind weights into local scope
        NumEqVector flux(0.0);

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                auto upwindTerm = [phaseIdx, compIdx](const VolumeVariables& volVars)
                { return volVars.molarDensity(phaseIdx)*volVars.moleFraction(phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                // get equation index
                auto eqIdx = Indices::conti0EqIdx + compIdx;
                flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
            }

            // Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        // Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        // diffusive fluxes
        const auto diffusionFluxesWPhase = fluxVars.molecularDiffusionFlux(wPhaseIdx);
        Scalar jGW = diffusionFluxesWPhase[gCompIdx];
        Scalar jNW = diffusionFluxesWPhase[nCompIdx];
        Scalar jWW = -(jGW+jNW);

        //check for the reference system and adapt units of the diffusive flux accordingly.
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
        {
            jGW /= FluidSystem::molarMass(gCompIdx);
            jNW /= FluidSystem::molarMass(nCompIdx);
            jWW /= FluidSystem::molarMass(wCompIdx);
        }

        const auto diffusionFluxesGPhase = fluxVars.molecularDiffusionFlux(gPhaseIdx);
        Scalar jWG = diffusionFluxesGPhase[wCompIdx];
        Scalar jNG = diffusionFluxesGPhase[nCompIdx];
        Scalar jGG = -(jWG+jNG);

        //check for the reference system and adapt units of the diffusive flux accordingly.
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
        {
            jWG /= FluidSystem::molarMass(wCompIdx);
            jNG /= FluidSystem::molarMass(nCompIdx);
            jGG /= FluidSystem::molarMass(gCompIdx);
        }

        // At the moment we do not consider diffusion in the NAPL phase
        const Scalar jWN = 0.0;
        const Scalar jGN = 0.0;
        const Scalar jNN = 0.0;

        flux[contiWEqIdx] += jWW+jWG+jWN;
        flux[contiNEqIdx] += jNW+jNG+jNN;
        flux[contiGEqIdx] += jGW+jGG+jGN;

        return flux;
    }
};

} // end namespace Dumux

#endif
