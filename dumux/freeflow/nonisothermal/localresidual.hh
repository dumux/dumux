// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeFlowEnergyLocalResidual
 */
#ifndef DUMUX_FREE_FLOW_ENERGY_LOCAL_RESIDUAL_HH
#define DUMUX_FREE_FLOW_ENERGY_LOCAL_RESIDUAL_HH

#include <cmath>

#include <dumux/discretization/method.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template<class GridGeometry, class FluxVariables, class DiscretizationMethod, bool enableEneryBalance, bool isCompositional>
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
                                                typename GridGeometry::DiscretizationMethod,
                                                enableEneryBalance, isCompositional>;

/*!
 * \ingroup FreeflowNIModel
 * \brief Specialization for isothermal models, does nothing
 */
template<class GridGeometry, class FluxVariables, class DiscretizationMethod, bool isCompositional>
class FreeFlowEnergyLocalResidualImplementation<GridGeometry, FluxVariables, DiscretizationMethod, false, isCompositional>
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
                                                DiscretizationMethods::Staggered,
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
                                                                                 fvGeometry,
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
                                                DiscretizationMethods::Staggered,
                                                true, true>
    : public FreeFlowEnergyLocalResidualImplementation<GridGeometry,
                                                       FluxVariables,
                                                       DiscretizationMethods::Staggered,
                                                       true, false>
{
    using ParentType = FreeFlowEnergyLocalResidualImplementation<GridGeometry,
                                                                 FluxVariables,
                                                                 DiscretizationMethods::Staggered,
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
            // define the upstream direction according to the sign of the diffusive flux
            using std::signbit;
            const bool insideIsUpstream = !signbit(diffusiveFlux[compIdx]);
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
