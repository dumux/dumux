// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup FacetCoupling
  * \brief Modified upwind scheme for models using the box scheme
  *        with coupling across element facets.
  */
#ifndef DUMUX_MULTIDOMAIN_FACET_BOX_UPWINDSCHEME_HH
#define DUMUX_MULTIDOMAIN_FACET_BOX_UPWINDSCHEME_HH

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
class BoxFacetCouplingUpwindScheme
{
public:
    // applies a simple upwind scheme to the precalculated advective flux
    template<class FluxVariables, class UpwindTermFunction, class Scalar>
    static Scalar apply(const FluxVariables& fluxVars,
                        const UpwindTermFunction& upwindTerm,
                        Scalar flux, int phaseIdx)
    {
        // TODO: pass this from outside?
        static const Scalar upwindWeight = getParam<Scalar>("Flux.UpwindWeight");

        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& scvf = fluxVars.scvFace();
        const auto& insideScv = fluxVars.fvGeometry().scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // check if this is an interior boundary
        if (scvf.interiorBoundary())
        {
            const auto& cm = fluxVars.problem().couplingManager();
            const auto& outsideVolVars = cm.getLowDimVolVars(fluxVars.element(), scvf);
            if (std::signbit(flux)) // if sign of flux is negative
                return flux*(upwindWeight*upwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
            else
                return flux*(upwindWeight*upwindTerm(insideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
        }
        else
        {
            const auto& outsideScv = fluxVars.fvGeometry().scv(scvf.outsideScvIdx());
            const auto& outsideVolVars = elemVolVars[outsideScv];
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
