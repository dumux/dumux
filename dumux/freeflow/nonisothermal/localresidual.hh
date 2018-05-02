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
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeFlowEnergyLocalResidual
 */
#ifndef DUMUX_FREE_FLOW_ENERGY_LOCAL_RESIDUAL_HH
#define DUMUX_FREE_FLOW_ENERGY_LOCAL_RESIDUAL_HH

#include <dumux/discretization/methods.hh>

namespace Dumux {

// forward declaration
template<class FVGridGeometry, class FluxVariables, DiscretizationMethod discMethod, bool enableEneryBalance>
class FreeFlowEnergyLocalResidualImplementation;

/*!
 * \ingroup FreeflowNIModel
 * \brief Element-wise calculation of the local residual for non-isothermal
 *        free-flow models
 */
template<class FVGridGeometry, class FluxVariables, bool enableEneryBalance>
using FreeFlowEnergyLocalResidual =
      FreeFlowEnergyLocalResidualImplementation<FVGridGeometry,
                                                FluxVariables,
                                                FVGridGeometry::discMethod,
                                                enableEneryBalance>;

/*!
 * \ingroup FreeflowNIModel
 * \brief Specialization for isothermal models, does nothing
 */
template<class FVGridGeometry, class FluxVariables, DiscretizationMethod discMethod>
class FreeFlowEnergyLocalResidualImplementation<FVGridGeometry, FluxVariables, discMethod, false>
{
public:

    //! do nothing for the isothermal case
    template <typename... Args>
    static void fluidPhaseStorage(Args&&... args)
    {}

    //! do nothing for the isothermal case
    template <typename... Args>
    static void heatConvectionFlux(Args&&... args)
    {}

    //! do nothing for the isothermal case
    template <typename... Args>
    static void heatFlux(Args&&... args)
    {}
};

/*!
 * \ingroup FreeflowNIModel
 * \brief Specialization for staggered non-isothermal models
 */
template<class FVGridGeometry, class FluxVariables>
class FreeFlowEnergyLocalResidualImplementation<FVGridGeometry,
                                                FluxVariables,
                                                DiscretizationMethod::staggered,
                                                true>
{
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using HeatConductionType = typename FluxVariables::HeatConductionType;

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
        using Indices = typename ElementVolumeVariables::VolumeVariables::Indices;

        bool isOutflow = false;
        if(scvf.boundary())
        {
            const auto bcTypes = problem.boundaryTypes(element, scvf);
            if(bcTypes.isOutflow(Indices::energyBalanceIdx))
                isOutflow = true;
        }

        auto upwindTerm = [](const auto& volVars) { return volVars.density() * volVars.enthalpy(); };
        flux[localEnergyBalanceIdx] += FluxVariables::advectiveFluxForCellCenter(elemVolVars,
                                                                            elemFaceVars,
                                                                            scvf,
                                                                            upwindTerm,
                                                                            isOutflow);
        if(!isOutflow)
            flux[localEnergyBalanceIdx] += HeatConductionType::flux(element,
                                                               fvGeometry,
                                                               elemVolVars,
                                                               scvf);
    }
};

} // end namespace Dumux

#endif
