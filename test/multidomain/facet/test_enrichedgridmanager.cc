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
 * \brief Tests the grid creator
 */
#include <config.h>
#include <iostream>

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/unused.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/common/parameters.hh>
#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/enrichedgridmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>
#include <dumux/multidomain/facet/vertexmapper.hh>

int main(int argc, char** argv) try
{
    // initialize MPI, finalize is done automatically on exit
    DUNE_UNUSED const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Dumux::Parameters::init(argc, argv);

    using BulkGrid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>;
    using EnrichedGrid = Dune::UGGrid<2>;
    using FacetGrid = Dune::FoamGrid<1, 2>;

    using BulkGridView = typename BulkGrid::LeafGridView;
    using EnrichedGridView = typename EnrichedGrid::LeafGridView;

    // Create both grid for the matrix and the facet domain
    using GridManager = Dumux::FacetCouplingGridManager<BulkGrid, FacetGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& facetGridView = gridManager.template grid<1>().leafGridView();
    const auto gridData = gridManager.getGridData()->template getSubDomainGridData<0>();

    using GridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings>;
    GridAdapter gridAdapter(gridManager.getEmbeddings());

    using VertexMapper = Dumux::EnrichedVertexDofMapper<BulkGridView>;
    VertexMapper vertexMapper(bulkGridView);
    vertexMapper.enrich(facetGridView, gridAdapter, true);

    // create enriched bulk grid
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>;
    ElementMapper elementMapper(bulkGridView, Dune::mcmgElementLayout());

    Dumux::EnrichedGridManager<EnrichedGrid> enrichedGridManager;
    enrichedGridManager.init(bulkGridView, elementMapper, vertexMapper, gridData);

    auto& enrichedGrid = enrichedGridManager.grid();
    const auto& enrichedGridView = enrichedGrid.leafGridView();
    const auto enrichedGridData = enrichedGridManager.getGridData();

    using EnrichedElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<EnrichedGridView>;
    EnrichedElementMapper enrichedElementMapper(enrichedGridView, Dune::mcmgElementLayout());

    Dune::VTKWriter<BulkGridView> bulkVtkWriter(bulkGridView);
    Dune::VTKWriter<EnrichedGridView> enrichedVtkWriter(enrichedGridView);

    // add element markers to output
    std::vector<int> bulkMarkers(bulkGridView.size(0));
    std::vector<int> enrichedMarkers(enrichedGridView.size(0));

    for (const auto& e : elements(bulkGridView))
        bulkMarkers[elementMapper.index(e)] = gridData->getElementDomainMarker(e);
    for (const auto& e : elements(enrichedGridView))
        enrichedMarkers[enrichedElementMapper.index(e)] = enrichedGridData->getElementDomainMarker(e);

    bulkVtkWriter.addCellData(bulkMarkers, "marker");
    enrichedVtkWriter.addCellData(enrichedMarkers, "marker");

    // add boundary markers to boundary elements
    std::vector<int> bulkBoundaryMarkers(bulkGridView.size(0));
    std::vector<int> enrichedBoundaryMarkers(enrichedGridView.size(0));

    for (const auto& e : elements(bulkGridView))
        for (const auto& is : intersections(bulkGridView, e))
            if (is.boundary() && gridData->wasInserted(is))
                bulkBoundaryMarkers[elementMapper.index(e)] = gridData->getBoundaryDomainMarker(is);
    for (const auto& e : elements(enrichedGridView))
        for (const auto& is : intersections(enrichedGridView, e))
            if (is.boundary() && enrichedGridData->wasInserted(is))
                enrichedBoundaryMarkers[enrichedElementMapper.index(e)] = enrichedGridData->getBoundaryDomainMarker(is);

    bulkVtkWriter.addCellData(bulkBoundaryMarkers, "boundaryMarker");
    enrichedVtkWriter.addCellData(enrichedBoundaryMarkers, "boundaryMarker");

    bulkVtkWriter.write("bulk");
    enrichedVtkWriter.write("bulk_enriched");

    std::cout << "Finished writing vtk files" << std::endl;

    return 0;

}
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
