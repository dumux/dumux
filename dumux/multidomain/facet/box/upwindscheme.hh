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
