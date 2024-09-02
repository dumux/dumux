// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Class representing the information context to assemble finite volume fluxes at a sub-control volume face
 */
#ifndef DUMUX_FLUX_FVFLUXCONTEXT_HH
#define DUMUX_FLUX_FVFLUXCONTEXT_HH

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Class representing the information context to assemble finite volume fluxes at a sub-control volume face
 *
 * \tparam Problem the problem type to solve (for boundary conditions)
 * \tparam FVElementGeometry the element geometry type
 * \tparam ElementVolumeVariables the element volume variables type
 * \tparam ElementFluxVariablesCache the element flux variables cache type
 */
template<class Problem,
         class FVElementGeometry,
         class ElementVolumeVariables,
         class ElementFluxVariablesCache>
class FVFluxContext
{
    using Element = typename FVElementGeometry::Element;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
public:

    FVFluxContext(const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const SubControlVolumeFace& scvf,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementFluxVariablesCache& elemFluxVarsCache)
    : problem_(problem)
    , element_(element)
    , fvGeometry_(fvGeometry)
    , scvf_(scvf)
    , elemVolVars_(elemVolVars)
    , elemFluxVarsCache_(elemFluxVarsCache)
    {}

    const Problem& problem() const
    { return problem_; }

    const Element& element() const
    { return element_; }

    const SubControlVolumeFace& scvf() const
    { return scvf_; }

    const FVElementGeometry& fvGeometry() const
    { return fvGeometry_; }

    const ElementVolumeVariables& elemVolVars() const
    { return elemVolVars_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return elemFluxVarsCache_; }

private:
    const Problem& problem_;
    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    const SubControlVolumeFace& scvf_;
    const ElementVolumeVariables& elemVolVars_;
    const ElementFluxVariablesCache& elemFluxVarsCache_;
};

} // end namespace Dumux

#endif
