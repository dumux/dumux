// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxFlux
 * \brief This file contains the data which is required to calculate
 *        dispersive fluxes.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_DISPERSION_FLUX_HH
#define DUMUX_DISCRETIZATION_BOX_DISPERSION_FLUX_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/traits.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/facetensoraverage.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod, ReferenceSystemFormulation referenceSystem>
class DispersionFluxImplementation;

/*!
 * \ingroup BoxFlux
 * \brief Specialization of a dispersion flux for the box method
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class DispersionFluxImplementation<TypeTag, DiscretizationMethods::Box, referenceSystem>
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

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numComponents = ModelTraits::numFluidComponents();

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
        ComponentFluxVector componentFlux(0.0);

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxVarsCache.shapeValues();

        // density interpolation
        Scalar rhoMassOrMole(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto rho = massOrMolarDensity(elemVolVars[scv], referenceSystem, phaseIdx);
            rhoMassOrMole += rho * shapeValues[scv.indexInElement()][0];
        }

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            // collect the dispersion tensor, the fluxVarsCache and the shape values
            const auto& dispersionTensor = [&]()
            {
                const auto& tensor =
                    ModelTraits::CompositionalDispersionModel::compositionalDispersionTensor(problem, scvf, fvGeometry,
                                                                                             elemVolVars, elemFluxVarsCache,
                                                                                             phaseIdx, compIdx);

                if (Dune::FloatCmp::eq(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor(), 1e-6))
                    return insideVolVars.extrusionFactor()*tensor;
                else
                {
                    const auto insideTensor = insideVolVars.extrusionFactor() * tensor;
                    const auto outsideTensor = outsideVolVars.extrusionFactor() * tensor;

                    // the resulting averaged dispersion tensor
                    return faceTensorAverage(insideTensor, outsideTensor, scvf.unitOuterNormal());
                }
            }();

            // the mole/mass fraction gradient
            Dune::FieldVector<Scalar, dimWorld> gradX(0.0);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto x = massOrMoleFraction(elemVolVars[scv], referenceSystem, phaseIdx, compIdx);
                gradX.axpy(x, fluxVarsCache.gradN(scv.indexInElement()));
            }

            // compute the dispersion flux
            componentFlux[compIdx] = -1.0 * rhoMassOrMole * vtmv(scvf.unitOuterNormal(), dispersionTensor, gradX)*Extrusion::area(fvGeometry, scvf);
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx))
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
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
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // collect the dispersion tensor
        const auto& dispersionTensor = [&]()
        {
            const auto& tensor =
                ModelTraits::ThermalDispersionModel::thermalDispersionTensor(problem, scvf, fvGeometry,
                                                                             elemVolVars, elemFluxVarsCache,
                                                                             phaseIdx);

            if (Dune::FloatCmp::eq(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor(), 1e-6))
                return insideVolVars.extrusionFactor()*tensor;
            else
            {
                const auto insideTensor = insideVolVars.extrusionFactor() * tensor;
                const auto outsideTensor = outsideVolVars.extrusionFactor() * tensor;

                // the resulting averaged dispersion tensor
                return faceTensorAverage(insideTensor, outsideTensor, scvf.unitOuterNormal());
            }
        }();

        // compute the temperature gradient with the shape functions
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, GridView::dimensionworld> gradTemp(0.0);
        for (auto&& scv : scvs(fvGeometry))
            gradTemp.axpy(elemVolVars[scv].temperature(), fluxVarsCache.gradN(scv.indexInElement()));

        // compute the heat conduction flux
        return -1.0*vtmv(scvf.unitOuterNormal(), dispersionTensor, gradTemp)*Extrusion::area(fvGeometry, scvf);
    }

};

} // end namespace Dumux

#endif
