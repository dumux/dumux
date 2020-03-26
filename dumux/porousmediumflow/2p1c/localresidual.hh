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
 * \ingroup TwoPOneCModel
 *
 * \copydoc Dumux::TwoPOneCLocalResidual
 */

#ifndef DUMUX_2P1C_LOCAL_RESIDUAL_HH
#define DUMUX_2P1C_LOCAL_RESIDUAL_HH

#include <dumux/porousmediumflow/immiscible/localresidual.hh>

namespace Dumux {
/*!
 * \ingroup TwoPOneCModel
 * \brief Element-wise calculation of the residual for the fully implicit two-phase one-component flow model.
 */
template<class TypeTag>
class TwoPOneCLocalResidual : public ImmiscibleLocalResidual<TypeTag>
{
    using ParentType = ImmiscibleLocalResidual<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    static const auto numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    //! Evaluate the storage term within a given scv
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);
        // Compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            storage[Indices::conti0EqIdx] +=
                volVars.porosity()
                * volVars.saturation(phaseIdx) * volVars.density(phaseIdx);

            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        // The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }

    //! Evaluate the fluxes over a face of a sub control volume
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        NumEqVector flux;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // The physical quantities for which we perform upwinding
            auto upwindTerm = [phaseIdx](const auto& volVars)
                              { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx); };

            flux[Indices::conti0EqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);

            // Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        // Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

        return flux;
    }
};

} // end namespace Dumux

#endif
