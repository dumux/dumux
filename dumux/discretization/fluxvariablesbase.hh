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

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief Base class for the flux variables
 *        Actual flux variables inherit from this class
 */
template<class TypeTag, class Implementation>
class FluxVariablesBase
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

    enum{ enableFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGlobalFluxVariablesCache) };
    enum{ isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    FluxVariablesBase() : problemPtr_(nullptr), scvFacePtr_(nullptr)
    {}

    void init(const Problem& problem,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace &scvFace)
    {
        problemPtr_ = &problem;
        elementPtr_ = &element;
        scvFacePtr_ = &scvFace;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;

        // update the stencil if needed
        if (!enableFluxVarsCache)
            stencil_ = asImp_().computeStencil(problem, element, fvGeometry, scvFace);
    }

    // when caching is enabled, get the stencil from the cache class
    template <typename T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), Stencil>::type& stencil() const
    { return problem().model().fluxVarsCache(scvFace()).stencil(); }

    // when caching is disabled, return the private stencil variable. The update(...) routine has to be called beforehand.
    template <typename T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), Stencil>::type& stencil()
    { return stencil_; }

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

    Stencil computeStencil(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    { DUNE_THROW(Dune::InvalidStateException, "computeStencil() routine is not provided by the implementation."); }

private:

    Implementation &asImp_()
    {
        assert(static_cast<Implementation*>(this) != 0);
        return *static_cast<Implementation*>(this);
    }

    const Implementation &asImp_() const
    {
        assert(static_cast<const Implementation*>(this) != 0);
        return *static_cast<const Implementation*>(this);
    }

    const Problem* problemPtr_;              //! Pointer to the problem
    const Element* elementPtr_;              //! Pointer to the element at hand
    const FVElementGeometry* fvGeometryPtr_;
    const SubControlVolumeFace* scvFacePtr_; //! Pointer to the sub control volume face for which the flux variables are created
    const ElementVolumeVariables* elemVolVarsPtr_;
    Stencil stencil_;                        //! The flux stencil
};
} // end namespace

#endif
