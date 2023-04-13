// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
template<class GridGeometry, class DiscretizationMethod>
class UpwindSchemeImpl;

/*!
 * \ingroup Flux
 * \brief The upwind scheme used for the advective fluxes.
 *        This depends on the chosen discretization method.
 */
template<class GridGeometry>
using UpwindScheme = UpwindSchemeImpl<GridGeometry, typename GridGeometry::DiscretizationMethod>;

namespace Detail {

//! returns the upwind factor which is multiplied to the advective flux across the given scvf
template<class ElemVolVars, class SubControlVolumeFace, class UpwindTermFunction, class Scalar>
Scalar upwindSchemeMultiplier(const ElemVolVars& elemVolVars,
                              const SubControlVolumeFace& scvf,
                              const UpwindTermFunction& upwindTerm,
                              Scalar flux, int phaseIdx)
{
    // TODO: pass this from outside?
    static const Scalar upwindWeight = getParamFromGroup<Scalar>(elemVolVars.gridVolVars().problem().paramGroup(), "Flux.UpwindWeight");

    const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
    const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

    using std::signbit;
    if (signbit(flux)) // if sign of flux is negative
        return (upwindWeight*upwindTerm(outsideVolVars)
                + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
    else
        return (upwindWeight*upwindTerm(insideVolVars)
                + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));
}

} // end namespace Detail

//! Upwind scheme for control-volume finite element methods (uses the standard scheme)
template<class GridGeometry, class DM>
class UpwindSchemeImpl<GridGeometry, DiscretizationMethods::CVFE<DM>>
{
public:
    //! applies a simple upwind scheme to the precalculated advective flux
    template<class FluxVariables, class UpwindTermFunction, class Scalar>
    static Scalar apply(const FluxVariables& fluxVars,
                        const UpwindTermFunction& upwindTerm,
                        Scalar flux, int phaseIdx)
    {
        return apply(fluxVars.elemVolVars(), fluxVars.scvFace(), upwindTerm, flux, phaseIdx);
    }

    //! applies a simple upwind scheme to the precalculated advective flux across the given scvf
    template<class ElemVolVars, class SubControlVolumeFace, class UpwindTermFunction, class Scalar>
    static Scalar apply(const ElemVolVars& elemVolVars,
                        const SubControlVolumeFace& scvf,
                        const UpwindTermFunction& upwindTerm,
                        Scalar flux, int phaseIdx)
    {
        return flux*multiplier(elemVolVars, scvf, upwindTerm, flux, phaseIdx);
    }

    //! returns the upwind factor which is multiplied to the advective flux across the given scvf
    template<class ElemVolVars, class SubControlVolumeFace, class UpwindTermFunction, class Scalar>
    static Scalar multiplier(const ElemVolVars& elemVolVars,
                             const SubControlVolumeFace& scvf,
                             const UpwindTermFunction& upwindTerm,
                             Scalar flux, int phaseIdx)
    {
        return Detail::upwindSchemeMultiplier(elemVolVars, scvf, upwindTerm, flux, phaseIdx);
    }
};

//! Upwind scheme for the cell-centered tpfa scheme
template<class GridGeometry>
class UpwindSchemeImpl<GridGeometry, DiscretizationMethods::CCTpfa>
{
    using GridView = typename GridGeometry::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    //! returns the upwind factor which is multiplied to the advective flux across the given scvf
    template<class ElemVolVars, class SubControlVolumeFace, class UpwindTermFunction, class Scalar>
    static Scalar multiplier(const ElemVolVars& elemVolVars,
                             const SubControlVolumeFace& scvf,
                             const UpwindTermFunction& upwindTerm,
                             Scalar flux, int phaseIdx)
    {
        static_assert(dim == dimWorld, "Multiplier cannot be computed on surface grids using cell-centered scheme!");
        return Detail::upwindSchemeMultiplier(elemVolVars, scvf, upwindTerm, flux, phaseIdx);
    }

    //! applies a simple upwind scheme to the precalculated advective flux across the given scvf
    template<class ElemVolVars, class SubControlVolumeFace, class UpwindTermFunction, class Scalar>
    static Scalar apply(const ElemVolVars& elemVolVars,
                        const SubControlVolumeFace& scvf,
                        const UpwindTermFunction& upwindTerm,
                        Scalar flux, int phaseIdx)
    {
        static_assert(dim == dimWorld, "This upwind scheme cannot be used on surface grids using cell-centered schemes!");
        return flux*multiplier(elemVolVars, scvf, upwindTerm, flux, phaseIdx);
    }

    // For surface and network grids (dim < dimWorld) we have to do a special upwind scheme
    template<class FluxVariables, class UpwindTermFunction, class Scalar>
    static Scalar apply(const FluxVariables& fluxVars,
                        const UpwindTermFunction& upwindTerm,
                        Scalar flux, int phaseIdx)
    {
        const auto& scvf = fluxVars.scvFace();
        const auto& elemVolVars = fluxVars.elemVolVars();

        // standard scheme
        if constexpr (dim == dimWorld)
            return apply(elemVolVars, scvf, upwindTerm, flux, phaseIdx);

        // on non-branching points the standard scheme works
        if (scvf.numOutsideScvs() == 1)
            return flux*Detail::upwindSchemeMultiplier(elemVolVars, scvf, upwindTerm, flux, phaseIdx);

        // handle branching points with a more complicated upwind scheme
        // we compute a flux-weighted average of all inflowing branches
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        Scalar branchingPointUpwindTerm = 0.0;
        Scalar sumUpwindFluxes = 0.0;

        // if the inside flux is positive (outflow) do fully upwind and return flux
        using std::signbit;
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
};

//! Upwind scheme for cell-centered mpfa schemes
template<class GridGeometry>
class UpwindSchemeImpl<GridGeometry, DiscretizationMethods::CCMpfa>
: public UpwindSchemeImpl<GridGeometry, DiscretizationMethods::CCTpfa> {};

} // end namespace Dumux

#endif
