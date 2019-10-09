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
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeFlowEnergyLocalResidual
 */
#ifndef DUMUX_FREE_FLOW_ENERGY_LOCAL_RESIDUAL_HH
#define DUMUX_FREE_FLOW_ENERGY_LOCAL_RESIDUAL_HH

#include <dumux/discretization/method.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template<class GridGeometry, class FluxVariables, DiscretizationMethod discMethod, bool enableEneryBalance, bool isCompositional>
class FreeFlowEnergyLocalResidualImplementation;

/*!
 * \ingroup FreeflowNIModel
 * \brief Element-wise calculation of the local residual for non-isothermal
 *        free-flow models
 */
template<class GridGeometry, class FluxVariables, bool enableEneryBalance, bool isCompositional>
using FreeFlowEnergyLocalResidual =
      FreeFlowEnergyLocalResidualImplementation<GridGeometry,
                                                FluxVariables,
                                                GridGeometry::discMethod,
                                                enableEneryBalance, isCompositional>;

/*!
 * \ingroup FreeflowNIModel
 * \brief Specialization for isothermal models, does nothing
 */
template<class GridGeometry, class FluxVariables, DiscretizationMethod discMethod, bool isCompositional>
class FreeFlowEnergyLocalResidualImplementation<GridGeometry, FluxVariables, discMethod, false, isCompositional>
{
public:

    //! do nothing for the isothermal case
    template <typename... Args>
    static void fluidPhaseStorage(Args&&... args)
    {}

    //! do nothing for the isothermal case
    template <typename... Args>
    static void heatFlux(Args&&... args)
    {}
};

/*!
 * \ingroup FreeflowNIModel
 * \brief Specialization for staggered one-phase, non-isothermal models
 */
template<class GridGeometry, class FluxVariables>
class FreeFlowEnergyLocalResidualImplementation<GridGeometry,
                                                FluxVariables,
                                                DiscretizationMethod::staggered,
                                                true, false>
{
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:

    //! The energy storage in the fluid phase
    template<class NumEqVector, class VolumeVariables>
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const VolumeVariables& volVars)
    {
        static constexpr auto localEnergyBalanceIdx = NumEqVector::dimension - 1;
        storage[localEnergyBalanceIdx] += volVars.density() * volVars.internalEnergy();
    }

    //! The convective and conductive heat fluxes in the fluid phase
    template<class NumEqVector, class Problem, class ElementVolumeVariables, class ElementFaceVariables>
    static void heatFlux(NumEqVector& flux,
                         const Problem& problem,
                         const Element &element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFaceVariables& elemFaceVars,
                         const SubControlVolumeFace& scvf)
    {
        static constexpr auto localEnergyBalanceIdx = NumEqVector::dimension - 1;

        auto upwindTerm = [](const auto& volVars) { return volVars.density() * volVars.enthalpy(); };
        flux[localEnergyBalanceIdx] += FluxVariables::advectiveFluxForCellCenter(problem,
                                                                                 elemVolVars,
                                                                                 elemFaceVars,
                                                                                 scvf,
                                                                                 upwindTerm);

        flux[localEnergyBalanceIdx] += FluxVariables::HeatConductionType::flux(problem,
                                                                               element,
                                                                               fvGeometry,
                                                                               elemVolVars,
                                                                               scvf);
    }
};

/*!
 * \ingroup FreeflowNIModel
 * \brief Specialization for staggered compositional, non-isothermal models
 */
template<class GridGeometry, class FluxVariables>
class FreeFlowEnergyLocalResidualImplementation<GridGeometry,
                                                FluxVariables,
                                                DiscretizationMethod::staggered,
                                                true, true>
    : public FreeFlowEnergyLocalResidualImplementation<GridGeometry,
                                                       FluxVariables,
                                                       DiscretizationMethod::staggered,
                                                       true, false>
{
    using ParentType = FreeFlowEnergyLocalResidualImplementation<GridGeometry,
                                                                 FluxVariables,
                                                                 DiscretizationMethod::staggered,
                                                                 true, false>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    //! The convective and conductive heat fluxes in the fluid phase
    template<class NumEqVector, class Problem, class ElementVolumeVariables, class ElementFaceVariables>
    static void heatFlux(NumEqVector& flux,
                         const Problem& problem,
                         const Element &element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFaceVariables& elemFaceVars,
                         const SubControlVolumeFace& scvf)
    {
        ParentType::heatFlux(flux, problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf);

        static constexpr auto localEnergyBalanceIdx = NumEqVector::dimension - 1;
        auto diffusiveFlux = FluxVariables::MolecularDiffusionType::flux(problem, element, fvGeometry, elemVolVars, scvf);
        for (int compIdx = 0; compIdx < FluxVariables::numComponents; ++compIdx)
        {
            const bool insideIsUpstream = scvf.directionSign() == sign(diffusiveFlux[compIdx]);
            const auto& upstreamVolVars = insideIsUpstream ? elemVolVars[scvf.insideScvIdx()] : elemVolVars[scvf.outsideScvIdx()];

            if (FluxVariables::MolecularDiffusionType::referenceSystemFormulation() == ReferenceSystemFormulation::massAveraged)
                flux[localEnergyBalanceIdx] += diffusiveFlux[compIdx] * upstreamVolVars.componentEnthalpy(compIdx);
            else
                flux[localEnergyBalanceIdx] += diffusiveFlux[compIdx] * upstreamVolVars.componentEnthalpy(compIdx)* elemVolVars[scvf.insideScvIdx()].molarMass(compIdx);

        }
    }
};

} // end namespace Dumux

#endif
