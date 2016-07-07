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
 * \brief Base class for the volume variables vector
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FLUXVARSCACHEVECTOR_HH
#define DUMUX_DISCRETIZATION_BOX_FLUXVARSCACHEVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag>
class BoxFluxVariablesCacheVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type update(Problem& problem)
    {
        problemPtr_ = &problem;
        fluxVarsCache_.resize(problem.model().fvGeometries().numScvf());
        for (const auto& element : elements(problem.gridView()))
        {
            // bind the geometries and volume variables to the element (all the elements in stencil)
            problem.model().fvGeometries_().bind(element);

            const auto& elementGeometry = element.geometry();
            const auto& localBasis = problem.model().fvGeometries().feLocalBasis(elementGeometry.type());
            const auto& fvGeometry = problem.model().fvGeometries(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                (*this)[scvf].update(problem, element, elementGeometry, localBasis, scvf);
            }
        }
    }

    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type update(Problem& problem)
    {
        problemPtr_ = &problem;
    }

    // Function is called by the BoxLocalJacobian prior to flux calculations on the element.
    // We assume the FVGeometries to be bound at this point
    void bind(const Element& element)
    {
        bindElement(element);
    }

    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type
    bindElement(const Element& element)
    {
        const auto& fvGeometry = problem_().model().fvGeometries(element);
        const auto& elementGeometry = element.geometry();
        const auto& localBasis = problem_().model().fvGeometries().feLocalBasis(elementGeometry.type());

        // temporary resizing of the cache
        const auto numScvf = fvGeometry.numScvf();
        fluxVarsCache_.resize(numScvf);
        IndexType localScvfIdx = 0;
        for (const auto& scvf : scvfs(fvGeometry))
            fluxVarsCache_[localScvfIdx++].update(problem_(), element, elementGeometry, localBasis, scvf);
    }

    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type
    bindElement(const Element& element) {}

    // access operators in the case of caching
    template <typename T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    // access operators in the case of no caching
    template <typename T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.indexInElement()]; }

    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.indexInElement()]; }

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;
};

} // end namespace

#endif
