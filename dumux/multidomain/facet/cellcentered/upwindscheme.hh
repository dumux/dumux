// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetCoupling
 * \brief Modified upwind scheme for models using cell-centered schemes
 *        with coupling across element facets.
 */
#ifndef DUMUX_MIXEDDIMENSION_FACET_CC_UPWINDSCHEME_HH
#define DUMUX_MIXEDDIMENSION_FACET_CC_UPWINDSCHEME_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief The upwind scheme used for the advective fluxes.
 *        This is a modified scheme for models involving coupling
 *        with a lower-dimensional domain across the element facets.
 */
template<class GridGeometry>
class CCFacetCouplingUpwindScheme
{
    using GridView = typename GridGeometry::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    // For surface and network grids (dim < dimWorld) we have to do a special upwind scheme
    template<class FluxVariables, class UpwindTermFunction, class Scalar, int d = dim, int dw = dimWorld>
    static typename std::enable_if<(d < dw), Scalar>::type
    apply(const FluxVariables& fluxVars,
          const UpwindTermFunction& upwindTerm,
          Scalar flux, int phaseIdx)
    {
        static const Scalar upwindWeight = getParam<Scalar>("Flux.UpwindWeight");

        // the volume variables of the inside sub-control volume
        const auto& scvf = fluxVars.scvFace();
        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // check if this is an interior boundary
        const auto& cm = fluxVars.problem().couplingManager();
        if (cm.isOnInteriorBoundary(fluxVars.element(), scvf))
        {
            const auto& outsideVolVars = cm.getLowDimVolVars(fluxVars.element(), scvf);
            if (std::signbit(flux))
                return flux*(upwindWeight*upwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
            else
                return flux*(upwindWeight*upwindTerm(insideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
        }
        else
        {
            // check if this is a branching point
            if (scvf.numOutsideScvs() > 1)
            {
                // more complicated upwind scheme
                // we compute a flux-weighted average of all inflowing branches
                Scalar branchingPointUpwindTerm = 0.0;
                Scalar sumUpwindFluxes = 0.0;

                // if the inside flux is positive (outflow) do fully upwind and return flux
                if (!std::signbit(flux))
                    return upwindTerm(insideVolVars)*flux;
                else
                    sumUpwindFluxes += flux;

                for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                {
                    // compute the outside flux
                    const auto& fvGeometry = fluxVars.fvGeometry();
                    const auto outsideScvIdx = scvf.outsideScvIdx(i);
                    const auto outsideElement = fvGeometry.gridGeometry().element(outsideScvIdx);
                    const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);

                    using AdvectionType = typename FluxVariables::AdvectionType;
                    const auto outsideFlux = AdvectionType::flux(fluxVars.problem(),
                                                                 outsideElement,
                                                                 fvGeometry,
                                                                 elemVolVars,
                                                                 flippedScvf,
                                                                 phaseIdx,
                                                                 fluxVars.elemFluxVarsCache());

                    if (!std::signbit(outsideFlux))
                        branchingPointUpwindTerm += upwindTerm(elemVolVars[outsideScvIdx])*outsideFlux;
                    else
                        sumUpwindFluxes += outsideFlux;
                }

                // the flux might be zero
                if (sumUpwindFluxes != 0.0)
                    branchingPointUpwindTerm /= -sumUpwindFluxes;
                else
                    branchingPointUpwindTerm = 0.0;

                // upwind scheme (always do fully upwind at branching points)
                // a weighting here would lead to an error since the derivation is based on a fully upwind scheme
                // TODO How to implement a weight of e.g. 0.5
                if (std::signbit(flux))
                    return flux*branchingPointUpwindTerm;
                else
                    return flux*upwindTerm(insideVolVars);
            }
            // non-branching points and boundaries
            else
            {
                const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
                if (std::signbit(flux))
                    return flux*(upwindWeight*upwindTerm(outsideVolVars)
                                 + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
                else
                    return flux*(upwindWeight*upwindTerm(insideVolVars)
                                 + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
            }
        }
    }

    // For grids with dim == dimWorld we use a simple upwinding scheme
    template<class FluxVariables, class UpwindTermFunction, class Scalar, int d = dim, int dw = dimWorld>
    static typename std::enable_if<(d == dw), Scalar>::type
    apply(const FluxVariables& fluxVars,
          const UpwindTermFunction& upwindTerm,
          Scalar flux, int phaseIdx)
    {
        static const Scalar upwindWeight = getParam<Scalar>("Flux.UpwindWeight");

        const auto& scvf = fluxVars.scvFace();
        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // check if this is an interior boundary
        const auto& cm = fluxVars.problem().couplingManager();
        if (cm.isOnInteriorBoundary(fluxVars.element(), scvf))
        {
            // upwind scheme
            const auto& outsideVolVars = cm.getLowDimVolVars(fluxVars.element(), scvf);
            if (std::signbit(flux))
                return flux*(upwindWeight*upwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
            else
                return flux*(upwindWeight*upwindTerm(insideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
        }
        else
        {
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            if (std::signbit(flux)) // if sign of flux is negative
                return flux*(upwindWeight*upwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
            else
                return flux*(upwindWeight*upwindTerm(insideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
        }
    }
};

} // end namespace Dumux

#endif
