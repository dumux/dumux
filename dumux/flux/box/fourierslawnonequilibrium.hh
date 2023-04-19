// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxFlux
 * \brief This file contains the data which is required to calculate
 *        energy fluxes due to molecular diffusion with Fourier's law.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_NONEQUILIBRIUM_HH
#define DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_NONEQUILIBRIUM_HH

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/flux/facetensoraverage.hh>

namespace Dumux {

// forward declaration
template <class TypeTag, class DiscretizationMethod>
class FouriersLawNonEquilibriumImplementation;

/*!
 * \ingroup BoxFlux
 * \brief Specialization of Fourier's Law for the box method for thermal nonequilibrium models.
 */
template <class TypeTag>
class FouriersLawNonEquilibriumImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using ThermalConductivityModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr auto numEnergyEqSolid = getPropValue<TypeTag, Properties::NumEnergyEqSolid>();
    static constexpr auto numEnergyEqFluid = getPropValue<TypeTag, Properties::NumEnergyEqFluid>();
    static constexpr auto numEnergyEq = numEnergyEqSolid + numEnergyEqFluid;
    static constexpr auto sPhaseIdx = ModelTraits::numFluidPhases();

public:
    /*!
     * \brief Returns the heat flux within a fluid or solid
     *        phase (in J/s) across the given sub-control volume face.
     */
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // get inside and outside diffusion tensors and calculate the harmonic mean
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];
        const auto computeLambda = [&](const auto& v){
            if constexpr (numEnergyEq == 1)
                return v.effectiveThermalConductivity();
            else if constexpr (numEnergyEqFluid == 1)
                return (phaseIdx != sPhaseIdx)
                        ? v.effectiveFluidThermalConductivity()
                        : v.effectiveSolidThermalConductivity();
            else
                return v.effectivePhaseThermalConductivity(phaseIdx);
        };

        auto insideLambda = computeLambda(insideVolVars);
        auto outsideLambda = computeLambda(outsideVolVars);

        // scale by extrusion factor
        insideLambda *= insideVolVars.extrusionFactor();
        outsideLambda *= outsideVolVars.extrusionFactor();

        // the resulting averaged diffusion tensor
        const auto lambda = faceTensorAverage(insideLambda, outsideLambda, scvf.unitOuterNormal());

        // evaluate gradTemp at integration point
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        Dune::FieldVector<Scalar, GridView::dimensionworld> gradTemp(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            // compute the temperature gradient with the shape functions
            if (phaseIdx < numEnergyEqFluid)
                gradTemp.axpy(elemVolVars[scv].temperatureFluid(phaseIdx), fluxVarsCache.gradN(scv.indexInElement()));
            else
               gradTemp.axpy(elemVolVars[scv].temperatureSolid(), fluxVarsCache.gradN(scv.indexInElement()));
        }

        // compute the heat conduction flux
        return -1.0*vtmv(scvf.unitOuterNormal(), lambda, gradTemp)*Extrusion::area(fvGeometry, scvf);
    }
};

} // end namespace Dumux

#endif
