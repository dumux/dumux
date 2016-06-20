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
#ifndef DUMUX_GEOMECHANICS_IMPLICIT_STRESSVARIABLESCACHE_HH
#define DUMUX_GEOMECHANICS_IMPLICIT_STRESSVARIABLESCACHE_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief The stress variables cache classes
 *        stores matrices to recover the discplacement on a scv face and stencil
 */
template<class TypeTag, typename DiscretizationMethod = void>
class StressVariablesCache
{};

// specialization for the Box Method
template<class TypeTag>
class StressVariablesCache<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::Box>::type >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MechanicalLawType = typename GET_PROP_TYPE(TypeTag, MechanicalLawType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvFace)
    {
        FluxVariables fluxVars;
        stencil_ = fluxVars.computeStencil(problem, scvFace);
    }

    const Stencil& stencil() const
    { return stencil_; }

private:
    Stencil stencil_;
};

// specialization for the cell centered tpfa method
template<class TypeTag>
class StressVariablesCache<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::CCTpfa>::type >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MechanicalLawType = typename GET_PROP_TYPE(TypeTag, MechanicalLawType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvFace)
    {
        FluxVariables fluxVars;
        stencil_ = fluxVars.computeStencil(problem, scvFace);
    }

    const Stencil& stencil() const
    { return stencil_; }

private:
    Stencil stencil_;
};

} // end namespace

#endif
