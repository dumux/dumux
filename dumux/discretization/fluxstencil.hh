// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief The flux stencil specialized for different discretization schemes
 */
#ifndef DUMUX_DISCRETIZATION_FLUXSTENCIL_HH
#define DUMUX_DISCRETIZATION_FLUXSTENCIL_HH

#include <vector>

#include <dune/common/reservedvector.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief The flux stencil specialized for different discretization schemes
 * \note There might be different stencils used for e.g. advection and diffusion for schemes
 *       where the stencil depends on variables. Also schemes might even have solution dependent
 *       stencil. However, we always reserve the stencil or all DOFs that are possibly involved
 *       since we use the flux stencil for matrix and assembly. This might lead to some zeros stored
 *       in the matrix.
 */
template<class FVElementGeometry, class DiscretizationMethod = typename FVElementGeometry::GridGeometry::DiscretizationMethod>
class FluxStencil;

/*
 * \ingroup Discretization
 * \brief Flux stencil specialization for the cell-centered tpfa scheme
 * \tparam FVElementGeometry The local view on the finite volume grid geometry
 */
template<class FVElementGeometry>
class FluxStencil<FVElementGeometry, DiscretizationMethods::CCTpfa>
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

public:
    //! Each cell I couples to a cell J always only via one face
    using ScvfStencilIForJ = Dune::ReservedVector<GridIndexType, 1>;

    //! The flux stencil type
    using Stencil = typename SubControlVolumeFace::Traits::GridIndexStorage;

    //! Returns the flux stencil
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
 * \tparam FVElementGeometry The local view on the finite volume grid geometry
 */
template<class FVElementGeometry>
class FluxStencil<FVElementGeometry, DiscretizationMethods::CCMpfa>
{
    using GridGeometry = typename FVElementGeometry::GridGeometry;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    // Use the stencil type of the primary interaction volume
    using NodalIndexSet = typename GridGeometry::GridIVIndexSets::DualGridIndexSet::NodalIndexSet;

public:
    //! We don't know yet how many faces couple to a neighboring element
    using ScvfStencilIForJ = std::vector<GridIndexType>;

    //! The flux stencil type
    using Stencil = typename NodalIndexSet::NodalGridStencilType;

    //! Returns the indices of the elements required for flux calculation on an scvf.
    static const Stencil& stencil(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf)
    {
        const auto& gridGeometry = fvGeometry.gridGeometry();

        // return the scv (element) indices in the interaction region
        if (gridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
            return gridGeometry.gridInteractionVolumeIndexSets().secondaryIndexSet(scvf).gridScvIndices();
        else
            return gridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(scvf).gridScvIndices();
    }
};

} // end namespace Dumux

#endif
