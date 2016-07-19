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
 * \brief Global flux variable cache
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GLOBAL_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_BOX_GLOBAL_FLUXVARSCACHE_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/box/elementfluxvariablescache.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag, bool EnableGlobalFluxVariablesCache>
class BoxGlobalFluxVariablesCache
{};

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag>
class BoxGlobalFluxVariablesCache<TypeTag, true>
{
    // the local class needs access to the problem
    friend BoxElementFluxVariablesCache<TypeTag, true>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

public:
    void update(Problem& problem)
    {
        problemPtr_ = &problem;
        fluxVarsCache_.resize(problem.gridView().size(0));
        for (const auto& element : elements(problem.gridView()))
        {
            auto eIdx = problem.elementMapper().index(element);
            // bind the geometries and volume variables to the element (all the elements in stencil)
            auto fvGeometry = localView(problem.model().globalFvGeometry());
            fvGeometry.bind(element);

            fluxVarsCache_[eIdx].resize(fvGeometry.numScvf());
            for (auto&& scvf : scvfs(fvGeometry))
            {
                this->get(eIdx, scvf.index()).update(problem, element, fvGeometry, scvf);
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const BoxGlobalFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    // access operator
    const FluxVariablesCache& get(IndexType eIdx, IndexType scvfIdx) const
    { return fluxVarsCache_[eIdx][scvfIdx]; }

    // access operator
    FluxVariablesCache& get(IndexType eIdx, IndexType scvfIdx)
    { return fluxVarsCache_[eIdx][scvfIdx]; }

    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag>
class BoxGlobalFluxVariablesCache<TypeTag, false>
{
    // the local class needs access to the problem
    friend BoxElementFluxVariablesCache<TypeTag, false>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

public:

    void update(Problem& problem)
    { problemPtr_ = &problem; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const BoxGlobalFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    // currently bound element
    const Problem* problemPtr_;
};

} // end namespace

#endif
