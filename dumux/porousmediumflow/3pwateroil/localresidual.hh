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
 * \ingroup ThreePWaterOilModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase three-component fully implicit model.
 */

#ifndef DUMUX_3P2CNI_LOCAL_RESIDUAL_HH
#define DUMUX_3P2CNI_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {
/*!
 * \ingroup ThreePWaterOilModel
 * \brief Element-wise calculation of the local residual for problems
 *        using the ThreePWaterOil fully implicit model.
 */
template<class TypeTag>
class ThreePWaterOilLocalResidual : public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
protected:
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

        conti0EqIdx = Indices::conti0EqIdx,//!< Index of the mass conservation equation for the water component
        conti1EqIdx = conti0EqIdx + 1,//!< Index of the mass conservation equation for the contaminant component

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,
    };

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param problem The problem
     *  \param scv The sub-control-volume
     *  \param volVars The volume variables
     */
     NumEqVector computeStorage(const Problem& problem,
                                const SubControlVolume& scv,
                                const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                    int eqIdx = (compIdx == wCompIdx) ? conti0EqIdx : conti1EqIdx;
                    storage[eqIdx] += volVars.porosity()
                                      * volVars.saturation(phaseIdx)
                                      * massOrMoleDensity(volVars, phaseIdx)
                                      * massOrMoleFraction(volVars, phaseIdx, compIdx);
            }

            //! The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        //! The energy storage in the solid matrix
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

        // get upwind weights into local scope
        NumEqVector flux(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                 const auto upwindTerm = [&massOrMoleDensity, &massOrMoleFraction, phaseIdx, compIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*massOrMoleFraction(volVars, phaseIdx, compIdx)*volVars.mobility(phaseIdx); };

                // get equation index
                auto eqIdx = conti0EqIdx + compIdx;
                flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
            }

            //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        // diffusive fluxes
        static constexpr auto referenceSystemFormulation = FluxVariables::MolecularDiffusionType::referenceSystemFormulation();
        const auto diffusionFluxesWPhase = fluxVars.molecularDiffusionFlux(wPhaseIdx);
        Scalar jNW = diffusionFluxesWPhase[nCompIdx];
        Scalar jWW = -jNW;
        // check for the reference system and adapt units of the diffusive flux accordingly.
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
        {
            jNW /=  FluidSystem::molarMass(nCompIdx);
            jWW /=  FluidSystem::molarMass(wCompIdx);
        }

        const auto diffusionFluxesGPhase = fluxVars.molecularDiffusionFlux(gPhaseIdx);
        Scalar jWG = diffusionFluxesGPhase[wCompIdx];
        Scalar jNG = -jWG;
        // check for the reference system and adapt units of the diffusive flux accordingly.
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
        {
            jWG /= FluidSystem::molarMass(wCompIdx);
            jNG /= FluidSystem::molarMass(nCompIdx);
        }

        const auto diffusionFluxesNPhase = fluxVars.molecularDiffusionFlux(nPhaseIdx);
        Scalar jWN = diffusionFluxesNPhase[wCompIdx];
        Scalar jNN = -jWN;
        // check for the reference system and adapt units of the diffusive flux accordingly.
        if (referenceSystemFormulation == ReferenceSystemFormulation::massAveraged)
        {
            jWN /= FluidSystem::molarMass(wCompIdx);
            jNN /= FluidSystem::molarMass(nCompIdx);
        }

        flux[conti0EqIdx] += jWW+jWG+jWN;
        flux[conti1EqIdx] += jNW+jNG+jNN;

        return flux;
    }
};

} // end namespace Dumux

#endif
