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
 * \brief Tests the grid manager class for models using facet coupling.
 */

#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/vtkfunction.hh>
#include <dumux/multidomain/facet/gridmanager.hh>

#ifndef BULKGRIDTYPE // default to ug grid if not provided by CMake
#define BULKGRIDTYPE Dune::UGGrid<3>
#endif

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse command line argument parameters
    Dumux::Parameters::init(argc, argv);

    // make grid from input file
    using BulkGrid = BULKGRIDTYPE;
    using FacetGrid = Dune::FoamGrid<2, 3>;
    using EdgeGrid = Dune::FoamGrid<1, 3>;
    using GridManager = Dumux::FacetCouplingGridManager<BulkGrid, FacetGrid, EdgeGrid>;
    GridManager gridManager;
    gridManager.init();

    // check grid sizes
    const auto& bulkGridView = gridManager.grid<0>().leafGridView();
    const auto& facetGridView = gridManager.grid<1>().leafGridView();
    const auto& edgeGridView = gridManager.grid<2>().leafGridView();

    if (bulkGridView.size(0) != 491)
        DUNE_THROW(Dune::InvalidStateException, "Bulk grid has " << bulkGridView.size(0) << " instead of 491 elements");
    if (bulkGridView.size(3) != 153)
        DUNE_THROW(Dune::InvalidStateException, "Bulk grid has " << bulkGridView.size(3) << " instead of 153 vertices");

    if (facetGridView.size(0) != 32)
        DUNE_THROW(Dune::InvalidStateException, "Facet grid has " << facetGridView.size(0) << " instead of 32 elements");
    if (facetGridView.size(2) != 23)
        DUNE_THROW(Dune::InvalidStateException, "Facet grid has " << facetGridView.size(2) << " instead of 23 vertices");

    if (edgeGridView.size(0) != 2)
        DUNE_THROW(Dune::InvalidStateException, "Edge grid has " << edgeGridView.size(0) << " instead of 2 elements");
    if (edgeGridView.size(1) != 3)
        DUNE_THROW(Dune::InvalidStateException, "Edge grid has " << edgeGridView.size(1) << " instead of 3 vertices");

    // all entities should have the domain index 10
    for (const auto& e : elements(bulkGridView))
    {
        const auto domainMarker = gridManager.getGridData()->getElementDomainMarker<0>(e);
        if (domainMarker != 10)
            DUNE_THROW(Dune::InvalidStateException, "Bulk element has domain marker " << domainMarker << " instead of 10");
    }

    for (const auto& e : elements(facetGridView))
    {
        const auto domainMarker = gridManager.getGridData()->getElementDomainMarker<1>(e);
        if (domainMarker != 10)
            DUNE_THROW(Dune::InvalidStateException, "Facet element has domain marker " << domainMarker << " instead of 10");
    }

    for (const auto& e : elements(edgeGridView))
    {
        const auto domainMarker = gridManager.getGridData()->getElementDomainMarker<2>(e);
        if (domainMarker != 10)
            DUNE_THROW(Dune::InvalidStateException, "Edge element has domain marker " << domainMarker << " instead of 10");
    }

    // Check boundary boundary segments. There should be mesh
    // file-defined segments on all sides except the positive x-direction
    std::size_t bSegmentCount = 0;
    std::size_t posXCount, negXCount, posYCount, negYCount, posZCount, negZCount;
    posXCount = negXCount = posYCount = negYCount = posZCount = negZCount = 0;
    for (const auto& e : elements(bulkGridView))
    {
        for (const auto& is : intersections(bulkGridView, e))
        {
            if (is.boundary())
            {
                const auto n = is.centerUnitOuterNormal();

                if (gridManager.getGridData()->wasInserted<0>(is))
                {
                    bSegmentCount++;
                    const auto insIdx = gridManager.getEmbeddings()->insertionIndex<0>(is);
                    const auto marker = gridManager.getGridData()->getBoundaryDomainMarker<0>(insIdx);

                    if (Dune::FloatCmp::eq(n[0], 1.0, 1e-6)) // pos x-dir
                    {
                        if (marker != 1)
                            DUNE_THROW(Dune::InvalidStateException, "positive x-face should have marker 7, but has " << marker);
                        posXCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[0], -1.0, 1e-6)) // neg x-dir
                    {
                        if (marker != 2)
                            DUNE_THROW(Dune::InvalidStateException, "negative x-face should have marker 4, but has " << marker);
                        negXCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[1], 1.0, 1e-6)) // pos y-dir
                    {
                        if (marker != 3)
                            DUNE_THROW(Dune::InvalidStateException, "positive y-face should have marker 5, but has " << marker);
                        posYCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[1], -1.0, 1e-6)) // neg y-dir
                    {
                        if (marker != 4)
                            DUNE_THROW(Dune::InvalidStateException, "negative y-face should have marker 6, but has " << marker);
                        negYCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[2], 1.0, 1e-6)) // pos z-dir
                    {
                        if (marker != 5)
                            DUNE_THROW(Dune::InvalidStateException, "positive z-face should have marker 1, but has " << marker);
                        posZCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[2], -1.0, 1e-6)) // neg z-dir
                    {
                        if (marker != 6)
                            DUNE_THROW(Dune::InvalidStateException, "negative z-face should have marker 2, but has " << marker);
                        negZCount++;
                    }
                    else
                        DUNE_THROW(Dune::InvalidStateException, "Could not deduce boundary segment orientation");
                }
                else
                    DUNE_THROW(Dune::InvalidStateException, "Boundary intersection was not inserted. Can't obtain boundary segment index");
            }
        }
    }

    //! make sure we found the right number of boundary segments
    if (bSegmentCount != 240) DUNE_THROW(Dune::InvalidStateException, "Found " << bSegmentCount << " instead of 200 boundary segments");
    if (posXCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << posXCount << " instead of 40 boundary segments in pos x-direction");
    if (negXCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << negXCount << " instead of 40 boundary segments in neg x-direction");
    if (posYCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << posYCount << " instead of 40 boundary segments in pos y-direction");
    if (negYCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << negYCount << " instead of 40 boundary segments in neg y-direction");
    if (posZCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << posZCount << " instead of 40 boundary segments in pos z-direction");
    if (negZCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << negZCount << " instead of 40 boundary segments in neg z-direction");

    //! check if we found the right number of embeddings
    const auto& edgeEmbeddings = gridManager.getEmbeddings()->embeddedEntityMap(2);
    const auto& edgeEmbedments = gridManager.getEmbeddings()->adjoinedEntityMap(2);
    if (edgeEmbeddings.size() != 0) DUNE_THROW(Dune::InvalidStateException, "The grid with lowest dimension can't have embedded entities");
    if (edgeEmbedments.size() != 2) DUNE_THROW(Dune::InvalidStateException, "Found " << edgeEmbedments.size() << " instead of 2 edge element embedments");
    for (const auto& embedments : edgeEmbedments)
        if (embedments.second.size() != 4)
            DUNE_THROW(Dune::InvalidStateException, "edge element is embedded in " << embedments.second.size() << " facet elements instead of 4");

    const auto& facetEmbeddings = gridManager.getEmbeddings()->embeddedEntityMap(1);
    const auto& facetEmbedments = gridManager.getEmbeddings()->adjoinedEntityMap(1);
    if (facetEmbeddings.size() != 8) DUNE_THROW(Dune::InvalidStateException, "Found " << facetEmbeddings.size() << " instead of 8 embeddings in facet grid");
    if (facetEmbedments.size() != 32) DUNE_THROW(Dune::InvalidStateException, "Found " << facetEmbedments.size() << " instead of 32 facet element embedments");
    for (const auto& embedments : facetEmbedments)
        if (embedments.second.size() != 2)
            DUNE_THROW(Dune::InvalidStateException, "facet element is embedded in " << embedments.second.size() << " bulk elements instead of 2");
    for (const auto& embeddings : facetEmbeddings)
        if (embeddings.second.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "facet element has " << embeddings.second.size() << " embedded entities instead of 1");

    const auto& bulkEmbeddings = gridManager.getEmbeddings()->embeddedEntityMap(0);
    const auto& bulkEmbedments = gridManager.getEmbeddings()->adjoinedEntityMap(0);
    if (bulkEmbeddings.size() != 56) DUNE_THROW(Dune::InvalidStateException, "Found " << bulkEmbeddings.size() << " instead of 56 embeddings in bulk grid");
    if (bulkEmbedments.size() != 0) DUNE_THROW(Dune::InvalidStateException, "The grid with highest dimension can't have embedments");

    std::size_t singleEmbeddings = 0;
    std::size_t doubleEmbeddings = 0;
    for (const auto& embeddings : bulkEmbeddings)
    {
        if (embeddings.second.size() != 1 && embeddings.second.size() != 2)
            DUNE_THROW(Dune::InvalidStateException, "bulk element has " << embeddings.second.size() << " embedded entities instead of 1 or 2");
        if (embeddings.second.size() == 1) singleEmbeddings++;
        if (embeddings.second.size() == 2) doubleEmbeddings++;
    }

    if (singleEmbeddings != 48) DUNE_THROW(Dune::InvalidStateException, "Found " << singleEmbeddings << " instead of 48 bulk elements with 1 embedding");
    if (doubleEmbeddings != 8) DUNE_THROW(Dune::InvalidStateException, "Found " << doubleEmbeddings << " instead of 8 bulk elements with 2 embeddings");

    // test access to embeddings for single elements
    singleEmbeddings = 0;
    doubleEmbeddings = 0;
    for (const auto& e : elements(bulkGridView))
    {
        if (gridManager.getEmbeddings()->embeddedEntityIndices<0>(e).size() == 1)
            singleEmbeddings++;
        else if (gridManager.getEmbeddings()->embeddedEntityIndices<0>(e).size() == 2)
            doubleEmbeddings++;
        else if (gridManager.getEmbeddings()->embeddedEntityIndices<0>(e).size() != 0)
            DUNE_THROW(Dune::InvalidStateException, "wrong number of embeddings!");

        if (gridManager.getEmbeddings()->adjoinedEntityIndices<0>(e).size() != 0)
            DUNE_THROW(Dune::InvalidStateException, "bulk grid can't be embedded anywhere!");
    }

    if (singleEmbeddings != 48) DUNE_THROW(Dune::InvalidStateException, "Found " << singleEmbeddings << " instead of 48 bulk elements with 1 embedding");
    if (doubleEmbeddings != 8) DUNE_THROW(Dune::InvalidStateException, "Found " << doubleEmbeddings << " instead of 8 bulk elements with 2 embeddings");

    std::size_t embeddings = 0;
    std::size_t embedments = 0;
    for (const auto& e : elements(facetGridView))
    {
        if (gridManager.getEmbeddings()->embeddedEntityIndices<1>(e).size() == 1)
            embeddings++;
        else if (gridManager.getEmbeddings()->embeddedEntityIndices<1>(e).size() != 0)
            DUNE_THROW(Dune::InvalidStateException, "wrong number of embeddings!");

        if (gridManager.getEmbeddings()->adjoinedEntityIndices<1>(e).size() == 2)
            embedments++;
        else if (gridManager.getEmbeddings()->adjoinedEntityIndices<1>(e).size() != 0)
            DUNE_THROW(Dune::InvalidStateException, "wrong number of embedments!");
    }

    if (embeddings != 8) DUNE_THROW(Dune::InvalidStateException, "Found " << embeddings << " instead of 8 embeddings in facet grid");
    if (embedments != 32) DUNE_THROW(Dune::InvalidStateException, "Found " << embedments << " instead of 32 facet element embedments");

    embeddings = 0;
    embedments = 0;
    for (const auto& e : elements(edgeGridView))
    {
        if (gridManager.getEmbeddings()->embeddedEntityIndices<2>(e).size() != 0)
            DUNE_THROW(Dune::InvalidStateException, "wrong number of embeddings!");

        if (gridManager.getEmbeddings()->adjoinedEntityIndices<2>(e).size() == 4)
            embedments++;
        else if (gridManager.getEmbeddings()->adjoinedEntityIndices<2>(e).size() != 0)
            DUNE_THROW(Dune::InvalidStateException, "wrong number of embedments!");
    }

    if (embeddings != 0) DUNE_THROW(Dune::InvalidStateException, "The grid with lowest dimension can't have embedded entities");
    if (embedments!= 2) DUNE_THROW(Dune::InvalidStateException, "Found " << edgeEmbedments.size() << " instead of 2 edge element embedments");

    // write .vtk file for each grid
    using BulkWriter = Dune::VTKWriter<typename BulkGrid::LeafGridView>;
    BulkWriter bulkWriter(bulkGridView); bulkWriter.write("bulkgrid");

    using FacetWriter = Dune::VTKWriter<typename FacetGrid::LeafGridView>;
    FacetWriter facetWriter(facetGridView); facetWriter.write("facetgrid");

    using EdgeWriter = Dune::VTKWriter<typename EdgeGrid::LeafGridView>;
    EdgeWriter edgeWriter(edgeGridView); edgeWriter.write("edgegrid");
}
