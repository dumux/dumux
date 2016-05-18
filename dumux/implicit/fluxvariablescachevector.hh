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
#ifndef DUMUX_IMPLICIT_FLUXVARSCACHEVECTOR_HH
#define DUMUX_IMPLICIT_FLUXVARSCACHEVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag>
class FluxVariablesCacheVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache)>::type update(const Problem& problem)
    {
        fluxVarsCache_.resize(problem.model().fvGeometries().numScvf());
        for (const auto& element : elements(problem.gridView()))
        {
            for (auto&& scvf : problem.model().fvGeometries(element).scvfs())
            {
                (*this)[scvf.index()].update(problem, element, scvf);
            }
        }
    }

    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableFluxVariablesCache)>::type update(const Problem& problem)
    {}

    template <typename T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache), FluxVariablesCache>::type& operator [](IndexType scvfIdx) const
    { return fluxVarsCache_[scvfIdx]; }

    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache), FluxVariablesCache>::type& operator [](IndexType scvfIdx)
    { return fluxVarsCache_[scvfIdx]; }

private:
    std::vector<FluxVariablesCache> fluxVarsCache_;
};

} // end namespace

#endif
