// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \brief Discharge through a sub-control-volume face under Manning's law
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_DISCHARGE_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_DISCHARGE_HH

#include <type_traits>

#include <dune/common/fvector.hh>

#include <dumux/discretization/extrusion.hh>

#include "approximation.hh"
#include "regularization.hh"

namespace Dumux::LongWave {

/*!
 * \ingroup LongWaveModel
 * \brief Discharge through a boundary face under normal flow, per unit area of the face
 *
 * The free surface is taken parallel to the bed, so the bed slope alone drives the flux and
 * the depth is the one of the inside sub-control volume. This describes an open outlet: water
 * leaves at the rate the local terrain conveys it, without backwater from downstream.
 */
template<class GradZ, class Normal, class Scalar>
Scalar normalFlowDischarge(const GradZ& gradZ, const Normal& unitOuterNormal,
                           const Scalar h, const Scalar manningN)
{ return -(gradZ*unitOuterNormal)*diffusivity(conveyance(h), manningN, gradZ.two_norm()); }

/*!
 * \ingroup LongWaveModel
 * \brief Volumetric discharge through a sub-control-volume face, positive when it leaves
 *        the inside sub-control volume of the face
 *
 * Manning's law driven by \f$ \nabla z + w \nabla h \f$, with the weight \f$ w \f$ selecting the
 * long-wave approximation (see `freeSurfaceWeight`), and the conveyance upwinded in flow direction.
 * The extrusion factor carries the width of a one-dimensional channel, so this is a volumetric
 * discharge in both one and two dimensions.
 *
 * The upwind weight can be set by the problem through `upwindWeight()` (default one, full
 * upwinding). A weight below one blends the conveyance of both sides, which removes the switch
 * of the derivative where the free-surface gradient across the face changes sign.
 *
 * This is the interior flux of the long-wave model. It is also what a gauge at an interior
 * cross-section measures, where no boundary condition can supply it.
 */
template<class Problem, class Element, class FVElementGeometry,
         class ElementVolumeVariables, class ElementFluxVariablesCache>
auto discharge(const Problem& problem,
               const Element& element,
               const FVElementGeometry& fvGeometry,
               const ElementVolumeVariables& elemVolVars,
               const typename FVElementGeometry::SubControlVolumeFace& scvf,
               const ElementFluxVariablesCache& elemFluxVarsCache)
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using Extrusion = Extrusion_t<GridGeometry>;
    using Scalar = std::decay_t<decltype(elemVolVars[scvf.insideScvIdx()].waterDepth())>;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    const auto n = problem.spatialParams().manningN(element);

    const auto& fluxVarCache = elemFluxVarsCache[scvf];
    const auto& shapeValues = fluxVarCache.shapeValues();
    Dune::FieldVector<Scalar, dimWorld> gradZ(0.0), gradh(0.0);
    Scalar h = 0.0;
    for (const auto& scv : scvs(fvGeometry))
    {
        const auto& volVars = elemVolVars[scv];
        const auto& gradN = fluxVarCache.gradN(scv.indexInElement());
        gradZ.axpy(volVars.bedSurface(), gradN);
        gradh.axpy(volVars.waterDepth(), gradN);
        h += volVars.waterDepth()*shapeValues[scv.indexInElement()][0];
    }

    auto gradH = gradZ;
    gradH.axpy(freeSurfaceWeight(gradZ.two_norm(), h, n), gradh);

    // all factors of the discharge except the free-surface gradient are positive, so the upwind
    // direction is known before the flux coefficient, which needs the upwind conveyance
    const auto gradHndA = -(gradH*scvf.unitOuterNormal())*Extrusion::area(fvGeometry, scvf)
                          *elemVolVars[scvf.insideScvIdx()].extrusionFactor();

    const auto& upwind = gradHndA > 0.0 ? elemVolVars[scvf.insideScvIdx()]
                                        : elemVolVars[scvf.outsideScvIdx()];
    const auto& downwind = gradHndA > 0.0 ? elemVolVars[scvf.outsideScvIdx()]
                                          : elemVolVars[scvf.insideScvIdx()];

    Scalar upwindWeight = 1.0;
    if constexpr (requires { problem.upwindWeight(); })
        upwindWeight = problem.upwindWeight();
    const auto hUpwind = upwindWeight*upwind.waterDepth() + (1.0 - upwindWeight)*downwind.waterDepth();

    return gradHndA*diffusivity(conveyance(hUpwind), n, gradH.two_norm());
}

} // end namespace Dumux::LongWave

#endif
