// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDFMModel
 * \brief Assembler for transfer fluxes across barriers in the box-dfm model.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_BOX_DFM_BARRIER_FLUXES_HH
#define DUMUX_POROUS_MEDIUM_FLOW_BOX_DFM_BARRIER_FLUXES_HH

#include <utility>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>

namespace Dumux {

//! Structure to hold barrier properties
template<class Scalar>
struct BoxDfmBarrierProperties
{
    Scalar aperture;
    Scalar normalPermeability;
};

/*!
 * \ingroup BoxDFMModel
 * \brief Assembler for transfer fluxes across barriers for immiscible porous medium flow models.
 */
template<class GridGeometry,
         class Problem,
         class FracturePropertyFunc>
class BoxDfmImmiscibleBarrierFluxes {
    using Scalar = typename Problem::Traits::Scalar;
    static_assert(
        int(GridGeometry::GridView::dimension) ==
        int(GridGeometry::GridView::dimensionworld)
    );

public:
    BoxDfmImmiscibleBarrierFluxes(const GridGeometry& gg,
                                  const Problem& p,
                                  FracturePropertyFunc&& fractureProperties)
    : gridGeometry_{gg}
    , problem_{p}
    , fracProperties_{std::move(fractureProperties)}
    {}

    template<class IntegrationPoint>
    auto operator()(const IntegrationPoint& ip) const {
        static_assert(
            std::is_invocable_v<FracturePropertyFunc, const IntegrationPoint&>,
            "Fracture property functor must be invocable with a barrier integration point"
        );
        static_assert(
            std::is_convertible_v<
                std::invoke_result_t<FracturePropertyFunc, const IntegrationPoint&>,
                BoxDfmBarrierProperties<Scalar>
            >,
            "Fracture property functor must return an instance of BoxDfmBarrierProperties<Scalar>"
        );
        const auto fractureProps = fracProperties_(ip);

        const auto& insideVolVars = ip.insideElemVolVars[ip.insideScv];
        const auto& outsideVolVars = ip.outsideElemVolVars[ip.outsideScv];
        const auto insideFVCache = makeFluxVarsCache_(ip.integrationPoint, ip.insideFVGeometry, ip.insideElemVolVars);
        const auto outsideFVCache = makeFluxVarsCache_(ip.integrationPoint, ip.outsideFVGeometry, ip.outsideElemVolVars);

        const auto& g = problem_.spatialParams().gravity(ip.integrationPoint);
        static const bool gravity = getParam<bool>("Problem.EnableGravity");

        typename Problem::Traits::NumEqVector result(0.0);
        for (int pIdx = 0; pIdx < result.size(); ++pIdx)
        {
            const auto insideP = getIPPressure_(pIdx, insideFVCache, ip.insideFVGeometry, ip.insideElemVolVars);
            const auto outsideP = getIPPressure_(pIdx, outsideFVCache, ip.outsideFVGeometry, ip.outsideElemVolVars);
            const auto gradP = (outsideP - insideP)/fractureProps.aperture;
            const auto rho = (insideVolVars.density(pIdx) + outsideVolVars.density(pIdx))/2.0;
            const auto flux = gravity
                ? -1.0*fractureProps.normalPermeability*(gradP - rho*(ip.normalVector*g))
                : -1.0*fractureProps.normalPermeability*gradP;
            result[pIdx] = flux > 0.0 ? flux*insideVolVars.mobility(pIdx) : flux*outsideVolVars.mobility(pIdx);
        }

        return result;
    }

private:
    template<class Pos, class FVGeometry, class ElementVolumeVariables>
    auto makeFluxVarsCache_(const Pos& ip, const FVGeometry& fvGeo, const ElementVolumeVariables& evv) const
    {
        Dumux::CVFEFluxVariablesCache<double, GridGeometry> cache;
        cache.update(problem_, fvGeo.element(), fvGeo, evv, ip);
        return cache;
    }

    template<class FVCache, class FVGeometry, class ElemVolVars>
    auto getIPPressure_(int phaseIdx,
                        const FVCache& cache,
                        const FVGeometry& fvGeo,
                        const ElemVolVars& evv) const
    {
        double result = 0.0;
        for (const auto& scv : scvs(fvGeo))
            if (!scv.isOnFracture())
                result += evv[scv].pressure(phaseIdx)*cache.shapeValues()[scv.localDofIndex()];
        return result;
    }

    const GridGeometry& gridGeometry_;
    const Problem& problem_;
    FracturePropertyFunc fracProperties_;
};

} // namespace Dumux

#endif
