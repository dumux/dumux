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
 * \brief Tests the grid creator class for models using facet coupling.
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
#include <dumux/multidomain/facet/gridcreator.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmapper.hh>

#ifndef BULKGRIDTYPE // default to ug grid if not provided by CMake
#define BULKGRIDTYPE Dune::UGGrid<3>
#endif

//! tests whether two positions are equal
template<typename Pos1, typename Pos2>
bool checkEquality(const Pos1& p1, const Pos2& p2, typename Pos1::value_type eps)
{
    const auto d = p1-p2;
    using std::abs;
    return std::all_of(d.begin(), d.end(), [eps] (auto coord) { return abs(coord) < eps; });
}

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse command line argument parameters
    Dumux::Parameters::init(argc, argv);

    using BulkGrid = BULKGRIDTYPE;
    using FacetGrid = Dune::FoamGrid<2, 3>;
    using EdgeGrid = Dune::FoamGrid<1, 3>;

    using GridCreator = Dumux::FacetCouplingGridCreator<BulkGrid, FacetGrid, EdgeGrid>;
    GridCreator gridCreator;
    gridCreator.makeGrids("grid.msh");

    // instantiate the grid geometries with caching
    using BulkGridView = typename BulkGrid::LeafGridView;
    using BulkFVGridGeometry = Dumux::CCTpfaFVGridGeometry<BulkGridView, true>;
    BulkFVGridGeometry bulkFvGeometry( gridCreator.grid<0>().leafGridView() );
    bulkFvGeometry.update();

    using FacetGridView = typename FacetGrid::LeafGridView;
    using FacetFVGridGeometry = Dumux::CCTpfaFVGridGeometry<FacetGridView, true>;
    FacetFVGridGeometry facetFvGeometry( gridCreator.grid<1>().leafGridView() );
    facetFvGeometry.update();

    using EdgeGridView = typename EdgeGrid::LeafGridView;
    using EdgeFVGridGeometry = Dumux::CCTpfaFVGridGeometry<EdgeGridView, true>;
    EdgeFVGridGeometry edgeFvGeometry( gridCreator.grid<2>().leafGridView() );
    edgeFvGeometry.update();

    // instantiate and update mappers for all domain combinations
    Dumux::CCTpfaFacetCouplingMapper<BulkFVGridGeometry, FacetFVGridGeometry> bulkFacetMapper;
    Dumux::CCTpfaFacetCouplingTwoDomainMapper<1, FacetFVGridGeometry, EdgeFVGridGeometry> facetEdgeMapper;
    Dumux::CCTpfaFacetCouplingMapper<BulkFVGridGeometry, FacetFVGridGeometry, EdgeFVGridGeometry> hierarchyMapper;

    bulkFacetMapper.update(bulkFvGeometry, facetFvGeometry, gridCreator);
    facetEdgeMapper.update(facetFvGeometry, edgeFvGeometry, gridCreator);
    hierarchyMapper.update(bulkFvGeometry, facetFvGeometry, edgeFvGeometry, gridCreator);

    constexpr auto bulkDomainId = Dune::index_constant<0>();
    constexpr auto facetDomainId = Dune::index_constant<1>();
    constexpr auto edgeDomainId = Dune::index_constant<2>();

    // check correctness of bulk-facet map
    for (unsigned int i = 0; i < 2; ++i)
    {
        // check both the map from the bulk facet as well as from the hierarchy mapper
        const auto& bulkFacetMap = i == 0 ? bulkFacetMapper.couplingMap(bulkDomainId, facetDomainId)
                                          : hierarchyMapper.couplingMap(bulkDomainId, facetDomainId);
        if (bulkFacetMap.size() != 56)
            DUNE_THROW(Dune::InvalidStateException, "BulkFacetMap has " << bulkFacetMap.size() << " instead of 56 entries");
        else
            std::cout << "Found 56 entries in bulk-facet map" << std::endl;

        std::size_t singleCouplings = 0;
        std::size_t doubleCouplings = 0;
        for (const auto& entry : bulkFacetMap)
        {
            const auto cStencilSize = entry.second.couplingStencil.size();

            if (cStencilSize == 1) singleCouplings++;
            else if (cStencilSize == 2) doubleCouplings++;
            else DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 1 or 2");

            for (unsigned int i = 0; i < cStencilSize; ++i)
            {
                const auto lowDimIdx = entry.second.couplingStencil[i];
                const auto bulkScvfIdx = entry.second.couplingScvfs[i][0];
                const auto lowDimGeom = facetFvGeometry.element(lowDimIdx).geometry();
                const auto& bulkScvf = bulkFvGeometry.scvf(bulkScvfIdx);
                if (!checkEquality(lowDimGeom.center(), bulkScvf.center(), lowDimGeom.volume()*1e-8))
                    DUNE_THROW(Dune::InvalidStateException, "Scvf does not coincide with low dim element");
            }
        }

        if (singleCouplings != 48) DUNE_THROW(Dune::InvalidStateException, "Found " << singleCouplings << " instead of 48 bulk coupling entries with size 1");
        if (doubleCouplings != 8) DUNE_THROW(Dune::InvalidStateException, "Found " << doubleCouplings << " instead of 8 bulk coupling entries with size 2");
    }

    // check correctness of facet-bulk map
    for (unsigned int i = 0; i < 2; ++i)
    {
        const auto& facetBulkMap = i == 0 ? bulkFacetMapper.couplingMap(facetDomainId, bulkDomainId)
                                          : hierarchyMapper.couplingMap(facetDomainId, bulkDomainId);
        if (facetBulkMap.size() != 32)
            DUNE_THROW(Dune::InvalidStateException, "FacetBulkMap has " << facetBulkMap.size() << " instead of 32 entries");
        else
            std::cout << "Found 32 entries in facet-bulk map" << std::endl;

        for (const auto& entry : facetBulkMap)
        {
            const auto cStencilSize = entry.second.couplingStencil.size();

            if (cStencilSize != 2)
                DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 2");

            for (unsigned int i = 0; i < cStencilSize; ++i)
            {
                const auto bulkIdx = entry.second.couplingStencil[i];
                const auto& embedment = entry.second.embedments[i];
                const auto lowDimGeom = facetFvGeometry.element(entry.first).geometry();
                const auto& bulkScvf = bulkFvGeometry.scvf(embedment.second[0]);
                if (bulkScvf.insideScvIdx() != bulkIdx)
                    DUNE_THROW(Dune::InvalidStateException, "Scvf insideScvIdx() does not match with bulk idx from embedment");
                if (!checkEquality(lowDimGeom.center(), bulkScvf.center(), lowDimGeom.volume()*1e-8))
                    DUNE_THROW(Dune::InvalidStateException, "Scvf does not coincide with low dim element");
            }
        }
    }

    // check correctness of facet-edge map
    for (unsigned int i = 0; i < 2; ++i)
    {
        const auto& facetEdgeMap = i == 0 ? facetEdgeMapper.couplingMap(facetDomainId, edgeDomainId)
                                          : hierarchyMapper.couplingMap(facetDomainId, edgeDomainId);
        if (facetEdgeMap.size() != 8)
            DUNE_THROW(Dune::InvalidStateException, "FacetEdgeMap has " << facetEdgeMap.size() << " instead of 8 entries");
        else
            std::cout << "Found 8 entries in facet-edge map" << std::endl;

        for (const auto& entry : facetEdgeMap)
        {
            const auto cStencilSize = entry.second.couplingStencil.size();

            if (cStencilSize != 1)
                DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 1");

            for (unsigned int i = 0; i < cStencilSize; ++i)
            {
                const auto lowDimIdx = entry.second.couplingStencil[i];
                const auto bulkScvfIdx = entry.second.couplingScvfs[i][0];
                const auto lowDimGeom = edgeFvGeometry.element(lowDimIdx).geometry();
                const auto& bulkScvf = facetFvGeometry.scvf(bulkScvfIdx);
                if (!checkEquality(lowDimGeom.center(), bulkScvf.center(), lowDimGeom.volume()*1e-8))
                    DUNE_THROW(Dune::InvalidStateException, "Scvf does not coincide with low dim element");
            }
        }
    }

    // check correctness of edge-facet map
    for (unsigned int i = 0; i < 2; ++i)
    {
        const auto& edgeFacetMap = i == 0 ? facetEdgeMapper.couplingMap(edgeDomainId, facetDomainId)
                                          : hierarchyMapper.couplingMap(edgeDomainId, facetDomainId);
        if (edgeFacetMap.size() != 2)
            DUNE_THROW(Dune::InvalidStateException, "EdgeFacetMap has " << edgeFacetMap.size() << " instead of 2 entries");
        else
            std::cout << "Found 2 entries in edge-facet map" << std::endl;

        for (const auto& entry : edgeFacetMap)
        {
            const auto cStencilSize = entry.second.couplingStencil.size();

            if (cStencilSize != 4)
                DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 4");

            for (unsigned int i = 0; i < cStencilSize; ++i)
            {
                const auto bulkIdx = entry.second.couplingStencil[i];
                const auto& embedment = entry.second.embedments[i];
                const auto lowDimGeom = edgeFvGeometry.element(entry.first).geometry();
                const auto& bulkScvf = facetFvGeometry.scvf(embedment.second[0]);
                if (bulkScvf.insideScvIdx() != bulkIdx)
                    DUNE_THROW(Dune::InvalidStateException, "Scvf insideScvIdx() does not match with bulk idx from embedment");
                if (!checkEquality(lowDimGeom.center(), bulkScvf.center(), lowDimGeom.volume()*1e-8))
                    DUNE_THROW(Dune::InvalidStateException, "Scvf does not coincide with low dim element");
            }
        }
    }

    // everything is ok
    std::cout << "\n... test passed!" << std::endl;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {
    std::cout << e << std::endl;
    return 1;
}
catch (...) {
    std::cout << "Unknown exception thrown" << std::endl;
    return 1;
}
