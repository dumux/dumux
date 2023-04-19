// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Base class for the flux variables living on a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH
#define DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH

#include <vector>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Base class for the flux variables living on a sub control volume face
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
class FluxVariablesBase
{
    using GridView = typename FVElementGeometry::GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Stencil = std::vector<std::size_t>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:

    //! Initialize the flux variables storing some temporary pointers
    void init(const Problem& problem,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace &scvFace,
              const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        problemPtr_ = &problem;
        elementPtr_ = &element;
        scvFacePtr_ = &scvFace;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        elemFluxVarsCachePtr_ = &elemFluxVarsCache;
    }

    const Problem& problem() const
    { return *problemPtr_; }

    const Element& element() const
    { return *elementPtr_; }

    const SubControlVolumeFace& scvFace() const
    { return *scvFacePtr_; }

    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return *elemFluxVarsCachePtr_; }

private:
    const Problem* problemPtr_;                             //!< Pointer to the problem
    const Element* elementPtr_;                             //!< Pointer to the element at hand
    const FVElementGeometry* fvGeometryPtr_;                //!< Pointer to the current FVElementGeometry
    const SubControlVolumeFace* scvFacePtr_;                //!< Pointer to the sub control volume face for which the flux variables are created
    const ElementVolumeVariables* elemVolVarsPtr_;          //!< Pointer to the current element volume variables
    const ElementFluxVariablesCache* elemFluxVarsCachePtr_; //!< Pointer to the current element flux variables cache
};

} // end namespace Dumux

#endif
