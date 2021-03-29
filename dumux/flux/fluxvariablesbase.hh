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
 * \brief Base class for the flux variables living on a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH
#define DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH

#include <vector>
#include <optional>

#include <dumux/discretization/method.hh>

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

namespace Experimental {

/*!
 * \ingroup Flux
 * \brief Base class for the flux variables living on a sub control volume face
 * \tparam LocalContext the element stencil-local context, consisting of
 *                      the local geometry and primary & secondary variables
 */
template<class LocalContext>
class FluxVariablesBase
{
    using FVElementGeometry = typename LocalContext::ElementGridGeometry;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = typename FVElementGeometry::GridGeometry;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Stencil = std::vector<std::size_t>;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;

public:

    //! Initialize the flux variables storing some temporary pointers
    void init(const LocalContext& context,
              const SubControlVolumeFace& scvf)
    {
        contextPtr_ = &context;
        scvFacePtr_ = &scvf;

        // for cell-centered methods, the element inside of the scvf may
        // not be the one the FVElementGeometry is bound to.
        if constexpr (!isBox)
        {
            const auto& fvGeometry = context.elementGridGeometry();
            const auto& gridGeometry = fvGeometry.gridGeometry();
            const auto& boundElement = fvGeometry.element();
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

            const auto boundElementIdx = gridGeometry.elementMapper().index(boundElement);
            if (insideScv.elementIndex() != boundElementIdx)
                element_ = gridGeometry.element(boundElementIdx);
        }
    }

    decltype(auto) problem() const
    { return elemVolVars().gridVolVars().problem(); }

    decltype(auto) elemVolVars() const
    { return contextPtr_->elementVariables().elemVolVars(); }

    decltype(auto) elemFluxVarsCache() const
    { return contextPtr_->elementVariables().elemFluxVarsCache(); }

    const SubControlVolumeFace& scvFace() const
    { return *scvFacePtr_; }

    const FVElementGeometry& fvGeometry() const
    { return contextPtr_->elementGridGeometry(); }

    const Element& element() const
    {
        if constexpr (isBox)
            return contextPtr_->elementGridGeometry().element();
        else
            return element_ ? *element_
                            : contextPtr_->elementGridGeometry().element();
    }

private:
    const LocalContext* contextPtr_;         //!< Pointer to the local context
    const SubControlVolumeFace* scvFacePtr_; //!< Pointer to the sub control volume face for which the flux variables are created
    std::optional<Element> element_;         //!< The element on the inside of the sub-control volume face
};

} // end namespace Experimental
} // end namespace Dumux

#endif
