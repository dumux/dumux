// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterModels
 * \copydoc Dumux::ShallowWaterResidual
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_LOCAL_RESIDUAL_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux{

/*!
 * \ingroup ShallowWaterModels
 * \brief Element-wise calculation of the residual for the shallow water equations
 */
template<class TypeTag>
class ShallowWaterResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;

public:

    using ParentType::ParentType;

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. mass, momentum) within a sub-control
     *        volume of a finite volume element.
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     * \note This function should not include the source and sink terms.
     * \note The volVars can be different to allow computing
     *       the implicit euler time derivative here
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        // partial time derivative of the phase mass
        NumEqVector storage(0.0);
        storage[Indices::massBalanceIdx] = volVars.waterDepth();
        storage[Indices::momentumXBalanceIdx] = volVars.waterDepth() * volVars.velocity(0);
        storage[Indices::momentumYBalanceIdx] = volVars.waterDepth() * volVars.velocity(1);
        return storage;
    }

       /*!
     * \brief Evaluate the mass/momentum flux over a face of a sub control volume
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache the flux variable cache for the element stencil
    */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        NumEqVector flux(0.0);
        FluxVariables fluxVars;
        flux += fluxVars.advectiveFlux(problem, element, fvGeometry, elemVolVars, scvf);

        // Compute viscous momentum flux contribution if required
        static const bool enableViscousFlux = getParamFromGroup<bool>(problem.paramGroup(), "ShallowWater.EnableViscousFlux", false);
        if (enableViscousFlux)
            flux += fluxVars.viscousFlux(problem, element, fvGeometry, elemVolVars, scvf);
        return flux;
    }
};
} // end namespace Dumux

#endif
