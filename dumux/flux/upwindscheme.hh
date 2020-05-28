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
 * \ingroup Flux
 * \brief Base class for the upwind scheme
 */
#ifndef DUMUX_DISCRETIZATION_UPWINDSCHEME_HH
#define DUMUX_DISCRETIZATION_UPWINDSCHEME_HH

#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

//! Forward declaration of the upwind scheme implementation
template<class GridGeometry, DiscretizationMethod discMethod>
class UpwindSchemeImpl;

/*!
 * \ingroup Flux
 * \brief The upwind scheme used for the advective fluxes.
 *        This depends on the chosen discretization method.
 */
template<class GridGeometry>
using UpwindScheme = UpwindSchemeImpl<GridGeometry, GridGeometry::discMethod>;

//! Upwind scheme for the box method
template<class GridGeometry>
class UpwindSchemeImpl<GridGeometry, DiscretizationMethod::box>
{
public:
    // applies a simple upwind scheme to the precalculated advective flux
    template<class FluxVariables, class UpwindTermFunction, class Scalar>
    static Scalar apply(const FluxVariables& fluxVars,
                        const UpwindTermFunction& upwindTerm,
                        Scalar flux, int phaseIdx)
    {
        // TODO: pass this from outside?
        static const Scalar upwindWeight = getParamFromGroup<Scalar>(fluxVars.problem().paramGroup(), "Flux.UpwindWeight");

        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& scvf = fluxVars.scvFace();
        const auto& insideScv = fluxVars.fvGeometry().scv(scvf.insideScvIdx());
        const auto& outsideScv = fluxVars.fvGeometry().scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        using std::signbit;
        if (signbit(flux)) // if sign of flux is negative
            return flux*(upwindWeight*upwindTerm(outsideVolVars)
                         + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
        else
            return flux*(upwindWeight*upwindTerm(insideVolVars)
                         + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
    }
};

//! Upwind scheme for the cell-centered tpfa scheme
template<class GridGeometry>
class UpwindSchemeImpl<GridGeometry, DiscretizationMethod::cctpfa>
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
        using std::signbit;
        static const Scalar upwindWeight = getParam<Scalar>("Flux.UpwindWeight");

        // the volume variables of the inside sub-control volume
        const auto& scvf = fluxVars.scvFace();
        const auto& elemVolVars = fluxVars.elemVolVars();
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        // check if this is a branching point
        if (scvf.numOutsideScvs() > 1)
        {
            // more complicated upwind scheme
            // we compute a flux-weighted average of all inflowing branches
            Scalar branchingPointUpwindTerm = 0.0;
            Scalar sumUpwindFluxes = 0.0;

            // if the inside flux is positive (outflow) do fully upwind and return flux
            if (!signbit(flux))
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

                if (!signbit(outsideFlux))
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
            if (signbit(flux))
                return flux*branchingPointUpwindTerm;
            else
                return flux*upwindTerm(insideVolVars);
        }
        // non-branching points and boundaries
        else
        {
            // upwind scheme
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
            if (signbit(flux))
                return flux*(upwindWeight*upwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
            else
                return flux*(upwindWeight*upwindTerm(insideVolVars)
                             + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
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
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        using std::signbit;
        if (signbit(flux)) // if sign of flux is negative
            return flux*(upwindWeight*upwindTerm(outsideVolVars)
                         + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
        else
            return flux*(upwindWeight*upwindTerm(insideVolVars)
                         + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
    }
};

//! Upwind scheme for cell-centered mpfa schemes
template<class GridGeometry>
class UpwindSchemeImpl<GridGeometry, DiscretizationMethod::ccmpfa>
: public UpwindSchemeImpl<GridGeometry, DiscretizationMethod::cctpfa> {};

} // end namespace Dumux

#endif
