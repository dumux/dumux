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
 * \ingroup Discretization
 * \brief Base class for the flux variables living on a sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH
#define DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH

#include <dumux/discretization/upwindscheme.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class Impl, class UpwindScheme>
class FluxVariablesBaseImplementation;

/*!
 * \ingroup Discretization
 * \brief Base class for the flux variables living on a sub control volume face
 * \note The upwind scheme is chosen depending on the discretization method
 */
template<class TypeTag, class Impl>
using FluxVariablesBase = FluxVariablesBaseImplementation<TypeTag, Impl, UpwindScheme<typename GET_PROP_TYPE(TypeTag, FVGridGeometry)>>;

/*!
 * \ingroup Discretization
 * \brief Base class for the flux variables living on a sub control volume face
 *
 * \tparam TypeTag The type tag
 * \tparam Implementation The implementation that uses this class as base (CRTP)
 * \tparam UpwindScheme The type of the upwind scheme used for upwinding of advective fluxes
 */
template<class TypeTag, class Implementation, class UpwindScheme>
class FluxVariablesBaseImplementation
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

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

    //! Applies the upwind scheme to precalculated fluxes
    template<class UpwindTermFunction>
    Scalar applyUpwindScheme(const UpwindTermFunction& upwindTerm, Scalar flux, int phaseIdx) const
    {
        //! Give the upwind scheme access to the cached variables
        return UpwindScheme::apply(asImp_(), upwindTerm, flux, phaseIdx);
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    const Problem* problemPtr_;                             //!< Pointer to the problem
    const Element* elementPtr_;                             //!< Pointer to the element at hand
    const FVElementGeometry* fvGeometryPtr_;                //!< Pointer to the current FVElementGeometry
    const SubControlVolumeFace* scvFacePtr_;                //!< Pointer to the sub control volume face for which the flux variables are created
    const ElementVolumeVariables* elemVolVarsPtr_;          //!< Pointer to the current element volume variables
    const ElementFluxVariablesCache* elemFluxVarsCachePtr_; //!< Pointer to the current element flux variables cache
};

} // end namespace Dumux

#endif
