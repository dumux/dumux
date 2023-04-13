// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \copydoc Dumux::ScheideggersDispersionTensor
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_SCHEIDEGGER_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_SCHEIDEGGER_HH

#include <algorithm>
#include <cmath>
#include <dune/common/math.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/facetensoraverage.hh>
#include <dumux/flux/traits.hh>

namespace Dumux {

namespace Detail {
template <class Problem, class SubControlVolumeFace>
using HasVelocityInSpatialParams = decltype(std::declval<Problem>().spatialParams().velocity(std::declval<SubControlVolumeFace>()));

template<class Problem, class SubControlVolumeFace>
static constexpr bool hasVelocityInSpatialParams()
{ return Dune::Std::is_detected<HasVelocityInSpatialParams, Problem, SubControlVolumeFace>::value; }
}

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Scheidegger's dispersion tensor
 */
template<class TypeTag>
class ScheideggersDispersionTensor
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static const int dimWorld = GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FluxTraits = typename Dumux::FluxTraits<FluxVariables>;
    static constexpr bool stationaryVelocityField = FluxTraits::hasStationaryVelocityField();

public:
    template <class ElementFluxVariablesCache>
    static DimWorldMatrix compositionalDispersionTensor(const Problem& problem,
                                                        const SubControlVolumeFace& scvf,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const ElementFluxVariablesCache& elemFluxVarsCache,
                                                        const int phaseIdx,
                                                        const int compIdx)
    {
        DimWorldMatrix dispersionTensor(0.0);

        // Get the velocity either from the reconstruction, or from the spatialparams
        auto velocity = dispersionVelocity_(problem, scvf, fvGeometry, elemVolVars, elemFluxVarsCache);

        // collect the dispersion alphas at this location
        std::array<Scalar,2> dispersivity = problem.spatialParams().dispersionAlphas(scvf.center(), phaseIdx, compIdx);

        return scheideggerTensor_(dispersivity, velocity);
    }

    template <class ElementFluxVariablesCache>
    static DimWorldMatrix thermalDispersionTensor(const Problem& problem,
                                                  const SubControlVolumeFace& scvf,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                                  const int phaseIdx)
    {
        DimWorldMatrix dispersionTensor(0.0);

        // Get the velocity either from the reconstruction, or from the spatialparams
        auto velocity = dispersionVelocity_(problem, scvf, fvGeometry, elemVolVars, elemFluxVarsCache);

        // collect the dispersion alphas at this location
        std::array<Scalar,2> dispersivity = problem.spatialParams().dispersionAlphas(scvf.center(), phaseIdx); //TODO: fix this?

        return scheideggerTensor_(dispersivity, velocity);
    }

private:

    template <class ElementFluxVariablesCache>
    static Dune::FieldVector<Scalar, dimWorld> dispersionVelocity_(const Problem& problem,
                                                                   const SubControlVolumeFace& scvf,
                                                                   [[maybe_unused]] const FVElementGeometry& fvGeometry,
                                                                   [[maybe_unused]] const ElementVolumeVariables& elemVolVars,
                                                                   [[maybe_unused]] const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // Calculate Darcy's velocity
        Dune::FieldVector<Scalar, dimWorld> velocity(0.0);
        if constexpr (stationaryVelocityField)
        {
            if constexpr (!Detail::hasVelocityInSpatialParams<Problem,SubControlVolumeFace>() )
                DUNE_THROW(Dune::NotImplemented, "\n Please provide the stationary velocity field in the spatialparams via a velocity function.");
            else
                velocity = problem.spatialParams().velocity(scvf);
        }
        else
        {
            if constexpr (FVElementGeometry::GridGeometry::discMethod == DiscretizationMethods::box)
            {
                const auto& fluxVarsCache = elemFluxVarsCache[scvf];
                const auto& shapeValues = fluxVarsCache.shapeValues();

                // get inside and outside permeability tensors and calculate the harmonic mean
                const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
                const auto K = faceTensorAverage(insideVolVars.permeability(),
                                                 outsideVolVars.permeability(),
                                                 scvf.unitOuterNormal());

                // evaluate gradP - rho*g at integration point
                Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
                Scalar rho(0.0);
                static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
                for (auto&& scv : scvs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[scv];

                    if (enableGravity)
                        rho += volVars.density(0)*shapeValues[scv.indexInElement()][0];
                    // the global shape function gradient
                    gradP.axpy(volVars.pressure(0), fluxVarsCache.gradN(scv.indexInElement()));
                }

                if (enableGravity)
                    gradP.axpy(-rho, problem.spatialParams().gravity(scvf.center()));

                velocity = gradP;
                velocity *= K;

                velocity /= -0.5 * (insideVolVars.viscosity() + outsideVolVars.viscosity());
            }
            else
                DUNE_THROW(Dune::NotImplemented, "\n Scheidegger Dispersion for compositional models without given constant velocity field is only implemented using the Box method.");
        }

        return velocity;
    }

    static DimWorldMatrix scheideggerTensor_(const std::array<Scalar,2>& dispersivity,
                                             const Dune::FieldVector<Scalar, dimWorld>& velocity)
    {
        DimWorldMatrix scheideggerTensor(0.0);

        //matrix multiplication of the velocity at the interface: vv^T
        for (int i=0; i < dimWorld; i++)
            for (int j = 0; j < dimWorld; j++)
                scheideggerTensor[i][j] = velocity[i]*velocity[j];

        //normalize velocity product --> vv^T/||v||, [m/s]
        Scalar vNorm = velocity.two_norm();

        scheideggerTensor /= vNorm;
        if (vNorm < 1e-20)
            scheideggerTensor = 0;

        //multiply with dispersivity difference: vv^T/||v||*(alphaL - alphaT), [m^2/s] --> alphaL = longitudinal disp., alphaT = transverse disp.
        scheideggerTensor *= (dispersivity[0] - dispersivity[1]);

        //add ||v||*alphaT to the main diagonal:vv^T/||v||*(alphaL - alphaT) + ||v||*alphaT, [m^2/s]
        for (int i = 0; i < dimWorld; i++)
            scheideggerTensor[i][i] += vNorm*dispersivity[1];

        return scheideggerTensor;
    }

};

} // end namespace Dumux

#endif
