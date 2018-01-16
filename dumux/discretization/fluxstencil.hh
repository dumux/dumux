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
 * \brief The flux stencil specialized for different discretization schemes
 */
#ifndef DUMUX_DISCRETIZATION_FLUXSTENCIL_HH
#define DUMUX_DISCRETIZATION_FLUXSTENCIL_HH

#include <dune/common/reservedvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

//! Forward declaration of the upwind scheme implementation
template<class TypeTag, DiscretizationMethods Method>
class FluxStencilImplementation;

/*!
 * \ingroup Discretization
 * \brief The flux stencil specialized for different discretization schemes
 * \note There might be different stencils used for e.g. advection and diffusion for schemes
 *       where the stencil depends on variables. Also schemes might even have solution dependent
 *       stencil. However, we always reserve the stencil or all DOFs that are possibly involved
 *       since we use the flux stencil for matrix and assembly. This might lead to some zeros stored
 *       in the matrix.
 */
template<class TypeTag>
using FluxStencil = FluxStencilImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

/*
 * \ingroup Discretization
 * \brief Flux stencil specialization for the cell-centered tpfa scheme
 */
template<class TypeTag>
class FluxStencilImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

public:
    //! The maximum number of elements in a flux stencil
    static constexpr int maxFluxStencilSize = GET_PROP_VALUE(TypeTag, MaxNumNeighborsPerScvf);

    //! Each cell I couples to a cell J always only via one face
    using ScvfStencilIForJ = Dune::ReservedVector<IndexType, 1>;

    //! The flux stencil type
    using Stencil = Dune::ReservedVector<IndexType, maxFluxStencilSize>;

    static Stencil stencil(const Element& element,
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

/*
 * \ingroup Discretization
 * \brief Flux stencil specialization for the cell-centered mpfa scheme
 */
template<class TypeTag>
class FluxStencilImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    static constexpr int dim = GridView::dimension;

    // Use the stencil type of the primary interaction volume
    using NodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);

public:
    //! The maximum number of elements in a flux stencil (equal to number of elements at node)
    static constexpr int maxFluxStencilSize = NodalIndexSet::maxNumElementsAtNode;

    //! We don't know yet how many faces couple to a neighboring element
    using ScvfStencilIForJ = std::vector<IndexType>;

    //! The flux stencil type
    using Stencil = typename NodalIndexSet::GridStencilType;

    /*
     * \brief Returns a set of grid element indices that participate in the
     *        flux calculations on a given scvf.
     *
     * \note The interaction volume index sets must use the same type for the
     *       stencils as the nodal index set. If not, the compiler will complain here.
     *
     * \param element The grid element
     * \param fvGeometry The finite volume geometry of this element
     * \param scvf The sub-control volume face embedded in this element
     */
    static const Stencil& stencil(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf)
    {
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();

        // return the scv (element) indices in the interaction region
        if (fvGridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
            return fvGridGeometry.gridInteractionVolumeIndexSets().secondaryIndexSet(scvf).globalScvIndices();
        else
            return fvGridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(scvf).globalScvIndices();
    }
};

} // end namespace Dumux

#endif
