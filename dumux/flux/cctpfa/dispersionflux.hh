// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCTpfaFlux
 * \brief This file contains the data which is required to calculate
 *        dispersive fluxes.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_DISPERSION_FLUX_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_DISPERSION_FLUX_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/traits.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod, ReferenceSystemFormulation referenceSystem>
class DispersionFluxImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Specialization of a Dispersion flux for the cctpfa method
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class DispersionFluxImplementation<TypeTag, DiscretizationMethods::CCTpfa, referenceSystem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridFluxVariablesCache = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVarCache = typename GridFluxVariablesCache::FluxVariablesCache;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FluxTraits = typename Dumux::FluxTraits<FluxVariables>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum
    {
        numPhases = ModelTraits::numFluidPhases(),
        numComponents = ModelTraits::numFluidComponents()
    };

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;
    using HeatFluxScalar = Scalar;

    static constexpr bool stationaryVelocityField = FluxTraits::hasStationaryVelocityField();

public:

    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    /*!
     * \brief Returns the dispersive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in mole/s or kg/s, depending
     *        on the template parameter ReferenceSystemFormulation.
     */
    static ComponentFluxVector compositionalDispersionFlux(const Problem& problem,
                                                           const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace& scvf,
                                                           const int phaseIdx,
                                                           const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (scvf.numOutsideScvs() > 1 )
            DUNE_THROW(Dune::NotImplemented, "\n Dispersion using ccTPFA is only implemented for conforming grids.");
        if (!stationaryVelocityField)
            DUNE_THROW(Dune::NotImplemented, "\n Dispersion using ccTPFA is only implemented for problems with stationary velocity fields");

        ComponentFluxVector componentFlux(0.0);
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto rhoInside = massOrMolarDensity(insideVolVars, referenceSystem, phaseIdx);
        const auto rhoOutside = massOrMolarDensity(outsideVolVars, referenceSystem, phaseIdx);
        const Scalar rho = 0.5*(rhoInside + rhoOutside);

        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            const auto& dispersionTensor =
                ModelTraits::CompositionalDispersionModel::compositionalDispersionTensor(problem, scvf, fvGeometry,
                                                                                         elemVolVars, elemFluxVarsCache,
                                                                                         phaseIdx, compIdx);
            const auto dij = computeTpfaTransmissibility(fvGeometry, scvf, fvGeometry.scv(scvf.insideScvIdx()), dispersionTensor, insideVolVars.extrusionFactor());

            const auto xInside = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const auto xOutide = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);

            componentFlux[compIdx] = (rho * (xInside-xOutide) * dij) * scvf.area();
        }
        return componentFlux;
    }

    /*!
     * \brief Returns the thermal dispersive flux
     *        across the given sub-control volume face.
     */
    static HeatFluxScalar thermalDispersionFlux(const Problem& problem,
                                                const Element& element,
                                                const FVElementGeometry& fvGeometry,
                                                const ElementVolumeVariables& elemVolVars,
                                                const SubControlVolumeFace& scvf,
                                                const int phaseIdx,
                                                const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (scvf.numOutsideScvs() > 1 )
            DUNE_THROW(Dune::NotImplemented, "\n Dispersion using ccTPFA is only implemented for conforming grids.");
        if (!stationaryVelocityField)
            DUNE_THROW(Dune::NotImplemented, "\n Dispersion using ccTPFA is only implemented for problems with stationary velocity fields");

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto& dispersionTensor =
            ModelTraits::ThermalDispersionModel::thermalDispersionTensor(problem, scvf, fvGeometry,
                                                                         elemVolVars, elemFluxVarsCache,
                                                                         phaseIdx);
        const auto dij = computeTpfaTransmissibility(scvf, fvGeometry.scv(scvf.insideScvIdx()), dispersionTensor, insideVolVars.extrusionFactor());

        // get the inside/outside temperatures
        const auto tInside = insideVolVars.temperature();
        const auto tOutside = outsideVolVars.temperature();

        // compute the heat conduction flux
        return (tInside-tOutside) * dij * scvf.area();
    }

};

} // end namespace Dumux

#endif
