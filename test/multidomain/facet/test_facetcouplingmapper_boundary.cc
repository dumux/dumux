// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup FacetTests
 * \brief Tests the coupling mapper class for models using facet coupling.
 */

#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dune/grid/uggrid.hh>
#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>

#ifndef BULKGRIDTYPE // default to ug grid if not provided by CMake
#define BULKGRIDTYPE Dune::UGGrid<2>
#endif

//! Tests whether two positions are equal.
template<typename Pos1, typename Pos2>
bool checkEquality(const Pos1& p1, const Pos2& p2, typename Pos1::value_type eps)
{
    const auto d = p1-p2;
    using std::abs;
    return std::all_of(d.begin(), d.end(), [eps] (auto coord) { return abs(coord) < eps; });
}

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Dumux::Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    using BulkGrid = BULKGRIDTYPE;
    using FacetGrid = Dune::FoamGrid<1, 2>;
    using GridManager = Dumux::FacetCouplingGridManager<BulkGrid, FacetGrid>;
    GridManager gridManager;
    gridManager.init();

    // instantiate the grid geometries with caching
    using BulkGridView = typename BulkGrid::LeafGridView;
    using BulkFVGridGeometry = Dumux::CCTpfaFVGridGeometry<BulkGridView, true>;
    BulkFVGridGeometry bulkFvGeometry( gridManager.grid<0>().leafGridView() );
    bulkFvGeometry.update();

    using FacetGridView = typename FacetGrid::LeafGridView;
    using FacetFVGridGeometry = Dumux::CCTpfaFVGridGeometry<FacetGridView, true>;
    FacetFVGridGeometry facetFvGeometry( gridManager.grid<1>().leafGridView() );
    facetFvGeometry.update();

    // instantiate and update mapper for all domain combinations
    Dumux::FacetCouplingMapper<BulkFVGridGeometry, FacetFVGridGeometry> mapper;
    mapper.update(bulkFvGeometry, facetFvGeometry, gridManager.getEmbeddings());

    constexpr auto bulkDomainId = Dune::index_constant<0>();
    constexpr auto facetDomainId = Dune::index_constant<1>();

    // check correctness of bulk-facet map
    const auto& bulkFacetMap = mapper.couplingMap(bulkDomainId, facetDomainId);
    if (bulkFacetMap.size() != 2)
        DUNE_THROW(Dune::InvalidStateException, "BulkFacetMap has " << bulkFacetMap.size() << " instead of 2 entries");
    else
        std::cout << "Found 2 entries in bulk-facet map" << std::endl;

    for (const auto& entry : bulkFacetMap)
    {
        const auto cStencilSize = entry.second.couplingStencil.size();

        if (cStencilSize != 1)
            DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 1");

        const auto lowDimIdx = entry.second.couplingStencil[0];
        const auto bulkScvfIdx = entry.second.dofToCouplingScvfMap.at(lowDimIdx)[0];
        const auto lowDimGeom = facetFvGeometry.element(lowDimIdx).geometry();
        const auto& bulkScvf = bulkFvGeometry.scvf(bulkScvfIdx);
        if (!checkEquality(lowDimGeom.center(), bulkScvf.center(), lowDimGeom.volume()*1e-8))
            DUNE_THROW(Dune::InvalidStateException, "Scvf does not coincide with low dim element");
    }

    // check correctness of facet-bulk map
    const auto& facetBulkMap = mapper.couplingMap(facetDomainId, bulkDomainId);
    if (facetBulkMap.size() != 2)
        DUNE_THROW(Dune::InvalidStateException, "FacetBulkMap has " << facetBulkMap.size() << " instead of 2 entries");
    else
        std::cout << "Found 2 entries in facet-bulk map" << std::endl;

    for (const auto& entry : facetBulkMap)
    {
        const auto cStencilSize = entry.second.couplingStencil.size();

        if (cStencilSize != 1)
            DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 1");

        const auto bulkIdx = entry.second.couplingStencil[0];
        const auto& embedment = entry.second.embedments[0];
        const auto lowDimGeom = facetFvGeometry.element(entry.first).geometry();
        const auto& bulkScvf = bulkFvGeometry.scvf(embedment.second[0]);
        if (bulkScvf.insideScvIdx() != bulkIdx)
            DUNE_THROW(Dune::InvalidStateException, "Scvf insideScvIdx() does not match with bulk idx from embedment");
        if (!checkEquality(lowDimGeom.center(), bulkScvf.center(), lowDimGeom.volume()*1e-8))
            DUNE_THROW(Dune::InvalidStateException, "Scvf does not coincide with low dim element");
    }

    // everything is ok
    std::cout << "\n... test passed!" << std::endl;
}
