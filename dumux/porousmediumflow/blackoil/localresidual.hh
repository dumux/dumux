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
 * \ingroup BlackOilModel
  * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the three-phase three-component fully implicit model.
 */

#ifndef DUMUX_BLACKOIL_LOCAL_RESIDUAL_HH
#define DUMUX_BLACKOIL_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux
{
/*!
 * \ingroup BlackOilModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the black-oil fully implicit model.
 */
template<class TypeTag>
class BlackOilLocalResidual : public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
/*
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
    */
    static constexpr int numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();
    static constexpr int numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();

   static constexpr int contiOEqIdx = Indices::conti0EqIdx + FluidSystem::wPhaseIdx;//!< index of the mass conservation equation for the water component
    static constexpr int contiWEqIdx = Indices::conti0EqIdx + FluidSystem::nPhaseIdx;//!< index of the mass conservation equation for the contaminant component
    static constexpr int contiGEqIdx = Indices::conti0EqIdx + FluidSystem::gPhaseIdx;//!< index of the mass conservation equation for the gas component

    static constexpr int wPhaseIdx = FluidSystem::wPhaseIdx;
    static constexpr int gPhaseIdx = FluidSystem::gPhaseIdx;
    static constexpr int oPhaseIdx = FluidSystem::oPhaseIdx;
    static constexpr int wCompIdx = FluidSystem::wCompIdx;
    static constexpr int oCompIdx = FluidSystem::oCompIdx; // NOTE: In the DUMUX setup the oily-phase is called "n"! In this model we changed this in the model!
    static constexpr int gCompIdx = FluidSystem::gCompIdx;

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
                                  * volVars.density(phaseIdx)
                                  * volVars.moleFraction(phaseIdx, compIdx);
            }
        }

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

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                auto upwindTerm = [phaseIdx, compIdx](const auto& volVars)
                {

                    return volVars.density(phaseIdx)*volVars.moleFraction(phaseIdx, compIdx)*volVars.mobility(phaseIdx);
                };

                // get equation index
                auto eqIdx = Indices::conti0EqIdx + compIdx;
                flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
            }

        }

        return flux;
    }
};

} // end namespace Dumux

#endif
