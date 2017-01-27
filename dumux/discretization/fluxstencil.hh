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
 * \brief Base class the flux stencil
 */
#ifndef DUMUX_DISCRETIZATION_FLUXSTENCIL_HH
#define DUMUX_DISCRETIZATION_FLUXSTENCIL_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

//! Forward declaration of the upwind scheme implementation
template<class TypeTag, DiscretizationMethods Method>
class FluxStencilImplementation;

/*!
 * \ingroup Discretization
 * \brief The flux stencil specialized for each discretization method
 * \note There might be different stencils used for e.g. advection and diffusion for schemes
 *       where the stencil depends on variables. Also schemes might even have solution dependent
 *       stencil. However, we always reserve the stencil or all DOFs that are possibly involved
 *       since we use the flux stencil for matrix and assembly. This might lead to some zeros stored
 *       in the matrix.
 */
template<class TypeTag>
using FluxStencil = FluxStencilImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

//! Flux stencil for the box method
template<class TypeTag>
class FluxStencilImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    // This is for compatibility with the cc methods. The flux stencil info is obsolete for the box method.
    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        return Stencil();
    }
};

//! Flux stencil for the cell-centered TPFA scheme
template<class TypeTag>
class FluxStencilImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        if (scvf.boundary())
            return Stencil({scvf.insideScvIdx()});
        else if (scvf.numOutsideScvs() > 1)
        {
            Stencil stencil({scvf.insideScvIdx()});
            for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                stencil.push_back(scvf.outsideScvIdx(i));
            return stencil;
        }
        else
            return Stencil({scvf.insideScvIdx(), scvf.outsideScvIdx()});
    }
};

//! Specialization for cell-centered MPFA schemes
template<class TypeTag>
class FluxStencilImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        const auto& globalFvGeometry = problem.model().globalFvGeometry();

        // return the scv (element) indices in the interaction region
        if (globalFvGeometry.touchesInteriorOrDomainBoundary(scvf))
            return globalFvGeometry.boundaryInteractionVolumeSeed(scvf).globalScvIndices();
        else
            return globalFvGeometry.interactionVolumeSeed(scvf).globalScvIndices();
    }
};

} // end namespace Dumux

#endif
