// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH
#define DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxstencil.hh>
#include <dumux/discretization/upwindscheme.hh>

namespace Dumux
{

template<class TypeTag, class UpwindScheme, class FluxStencil>
class FluxVariablesBaseImplementation;

/*!
 * \ingroup Discretization
 * \brief The flux variables base class class
 *        The upwind scheme is chosen depending on the discretization method
 */
template<class TypeTag>
using FluxVariablesBase = FluxVariablesBaseImplementation<TypeTag, UpwindScheme<TypeTag>, FluxStencil<TypeTag>>;

/*!
 * \ingroup Discretization
 * \brief Implementation of the base class of the flux variables
 *
 * \param TypeTag The type tag
 * \param UpwindScheme The type used for the upwinding of the advective fluxes
 */
template<class TypeTag, class UpwindScheme, class FluxStencil>
class FluxVariablesBaseImplementation
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

public:

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

    //! Applies the upwind scheme to precalculated fluxes
    template<class UpwindTermFunction>
    Scalar applyUpwindScheme(const UpwindTermFunction& upwindTerm, Scalar flux, int phaseIdx) const
    {
        //! Give the upwind scheme access to the cached variables
        return UpwindScheme::apply(*this, upwindTerm, flux, phaseIdx);
    }

    static Stencil computeStencil(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf)
    {
        //! Give the upwind scheme access to the cached variables
        //! Forward to the discretization specific implementation
        return FluxStencil::stencil(element, fvGeometry, scvf);
    }

private:
    const Problem* problemPtr_;              //! Pointer to the problem
    const Element* elementPtr_;              //! Pointer to the element at hand
    const FVElementGeometry* fvGeometryPtr_;
    const SubControlVolumeFace* scvFacePtr_; //! Pointer to the sub control volume face for which the flux variables are created
    const ElementVolumeVariables* elemVolVarsPtr_;
    const ElementFluxVariablesCache* elemFluxVarsCachePtr_;
};

} // end namespace Dumux

#endif
