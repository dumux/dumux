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
 * \ingroup Fluidmatrixinteractions
 * \copydoc Dumux::ScheideggersDispersionTensor
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_SCHEIDEGGER_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_DISPERSIONTENSORS_SCHEIDEGGER_HH

#include <algorithm>
#include <cmath>
#include <dune/common/math.hh>
#include <dumux/flux/facetensoraverage.hh>

namespace Dumux {

namespace Detail {
template <class Problem, class SubControlVolumeFace>
using HasVelocityInSpatialParams = decltype(std::declval<Problem>().spatialParams().velocity(std::declval<SubControlVolumeFace>()));

template<class Problem, class SubControlVolumeFace>
static constexpr bool hasVelocityInSpatialParams()
{ return Dune::Std::is_detected<HasVelocityInSpatialParams, Problem, SubControlVolumeFace>::value; }
}

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
                                                        [[maybe_unused]] const FVElementGeometry& fvGeometry,
                                                        [[maybe_unused]] const ElementVolumeVariables& elemVolVars,
                                                        [[maybe_unused]] const ElementFluxVariablesCache& elemFluxVarsCache,
                                                        [[maybe_unused]] const int phaseIdx,
                                                        [[maybe_unused]] const int compIdx)
    {
        DimWorldMatrix dispersionTensor(0.0);


    template <class ElementFluxVariablesCache>
    static DimWorldMatrix thermalDispersionTensor(const Problem& problem,
                                                  const SubControlVolumeFace& scvf,
                                                  [[maybe_unused]] const FVElementGeometry& fvGeometry,
                                                  [[maybe_unused]] const ElementVolumeVariables& elemVolVars,
                                                  [[maybe_unused]] const ElementFluxVariablesCache& elemFluxVarsCache,
                                                  [[maybe_unused]] const int phaseIdx)
    {
        DimWorldMatrix dispersionTensor(0.0);

        // Get the velocity either from the reconstruction, or from the spatialparams
        auto velocity = dispersionVelocity_(problem, scvf, fvGeometry, elemVolVars, elemFluxVarsCache);

        // collect the dispersion alphas at this location
        std::array<Scalar,2> dispersivity = problem.spatialParams().dispersionAlphas(scvf.center(), phaseIdx); //TODO: fix this?

        return scheideggerTensor_(dispersivity, velocity);
    }

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

        //calculate dispersion tensor

        // collect the dispersion alphas at this location
        std::array<Scalar,2> dispersivity = problem.spatialParams().dispersionAlphas(scvf.center(), phaseIdx, compIdx);

        //matrix multiplication of the velocity at the interface: vv^T
        for (int i=0; i < dimWorld; i++)
            for (int j = 0; j < dimWorld; j++)
                dispersionTensor[i][j] = velocity[i]*velocity[j];

        //normalize velocity product --> vv^T/||v||, [m/s]
        Scalar vNorm = velocity.two_norm();

        dispersionTensor /= vNorm;
        if (vNorm < 1e-20)
            dispersionTensor = 0;

        //multiply with dispersivity difference: vv^T/||v||*(alphaL - alphaT), [m^2/s] --> alphaL = longitudinal disp., alphaT = transverse disp.
        dispersionTensor *= (dispersivity[0] - dispersivity[1]);

        //add ||v||*alphaT to the main diagonal:vv^T/||v||*(alphaL - alphaT) + ||v||*alphaT, [m^2/s]
        for (int i = 0; i < dimWorld; i++)
            dispersionTensor[i][i] += vNorm*dispersivity[1];

        return dispersionTensor;
    }

};

} // end namespace Dumux

#endif
