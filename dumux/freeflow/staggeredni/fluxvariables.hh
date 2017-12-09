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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_FREELOW_IMPLICIT_NI_FLUXVARIABLES_HH
#define DUMUX_FREELOW_IMPLICIT_NI_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

namespace Properties
{
    NEW_PROP_TAG(ElementFaceVariables);
}

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables class
 *        specializations are provided for combinations of physical processes
 * \note  Not all specializations are currently implemented
 */

// forward declaration
template<class TypeTag, bool enableEnergyBalance>
class FreeFlowEnergyFluxVariablesImplementation;

template<class TypeTag>
using FreeFlowEnergyFluxVariables = FreeFlowEnergyFluxVariablesImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

// specialization for isothermal flow
template<class TypeTag>
class FreeFlowEnergyFluxVariablesImplementation<TypeTag, false>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);

public:

    static void energyFlux(CellCenterPrimaryVariables& flux,
                           const Problem& problem,
                           const Element &element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFaceVariables& elemFaceVars,
                           const SubControlVolumeFace &scvf,
                           const FluxVariablesCache& fluxVarsCache)
    { }

};

// specialization for non-isothermal flow
template<class TypeTag>
class FreeFlowEnergyFluxVariablesImplementation<TypeTag, true>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);

    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    enum { energyBalanceIdx = Indices::energyBalanceIdx };

public:

    static void energyFlux(CellCenterPrimaryVariables& flux,
                           const Problem& problem,
                           const Element &element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFaceVariables& elemFaceVars,
                           const SubControlVolumeFace &scvf,
                           const FluxVariablesCache& fluxVarsCache)
    {
        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        // TODO: catch neumann and outflow in localResidual's evalBoundary_()
        bool isOutflow = false;
        if(scvf.boundary())
        {
            const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
                if(bcTypes.isOutflow(energyBalanceIdx))
                    isOutflow = true;
        }

        auto upwindTerm = [](const auto& volVars) { return volVars.density() * volVars.enthalpy(); };

        flux[energyBalanceIdx] = FluxVariables::advectiveFluxForCellCenter(elemVolVars, elemFaceVars, scvf, upwindTerm, isOutflow);
        flux[energyBalanceIdx] += HeatConductionType::diffusiveFluxForCellCenter(problem, element, fvGeometry, elemVolVars, scvf);
    }

};

} // end namespace

#endif
