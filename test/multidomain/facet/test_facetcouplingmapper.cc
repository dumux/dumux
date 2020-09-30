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
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dune/grid/uggrid.hh>
#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometrytraits.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>
#include <dumux/multidomain/facet/box/fvgridgeometry.hh>
#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

#ifndef BULKGRIDTYPE // default to ug grid if not provided by CMake
#define BULKGRIDTYPE Dune::UGGrid<3>
#endif

#ifndef USEBOXINBULK // default to tpfa if not specified otherwise
#define USEBOXINBULK 0
#endif

#ifndef USEMPFAINBULK // default to tpfa if not specified otherwise
#define USEMPFAINBULK 0
#endif

//! Helper class to define the MpfaFVGridGeometry
template<class GridView>
struct MpfaGridGeometryTraits
{
private:
    using NITraits = Dumux::NodalIndexSetDefaultTraits< GridView >;
    using IndexSet = Dumux::CCMpfaDualGridNodalIndexSet< NITraits >;
    using IVTraits = Dumux::CCMpfaODefaultInteractionVolumeTraits< IndexSet, double >;
    using IV = Dumux::CCMpfaOInteractionVolume< IVTraits >;
    using GGTraits = Dumux::CCMpfaFVGridGeometryTraits<GridView, IndexSet, IV, IV>;
public:
    using type = Dumux::CCMpfaFVGridGeometry<GridView, GGTraits, /*enableCache*/true>;
};

//! Computes the average distance of corners to the center of a geometry
template<class Geometry>
typename Geometry::ctype averageCornerDistance(const Geometry& geometry)
{
    const auto center = geometry.center();
    typename Geometry::ctype avgDistance = 0.0;
    for (int i = 0; i < geometry.corners(); ++i)
        avgDistance += (geometry.corner(i) - center).two_norm();
    avgDistance /= double(geometry.corners());
    return avgDistance;
}

//! Checks if an scvf lies on a low dim element geometry.
template<class Scvf, class LowDimGeom>
void checkScvfEmbedment(const Scvf& scvf, const LowDimGeom& lowDimGeom)
{
    // scalar product of scvf center minus low dim element center and normal should be zero!
    const auto d = lowDimGeom.center()-scvf.center();
    const auto sp = d*scvf.unitOuterNormal();

    using std::abs;
    if ( !(abs(sp) < lowDimGeom.volume()*1e-15) || d.two_norm() > averageCornerDistance(lowDimGeom) )
        DUNE_THROW(Dune::InvalidStateException, "Scvf does not coincide with low dim element");
}

//! Updates a tpfa finite volume grid geometry.
template< class BulkFVG, class FacetFVG, class GridManager,
          std::enable_if_t<BulkFVG::discMethod != Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFvGeometry(BulkFVG& bulkFVG, const FacetFVG& facetFVG, const GridManager& gm)
{
    bulkFVG.update();
}

//! Updates a box finite volume grid geometry.
template< class BulkFVG, class FacetFVG, class GridManager,
          std::enable_if_t<BulkFVG::discMethod == Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFvGeometry(BulkFVG& bulkFVG, const FacetFVG& facetFVG, const GridManager& gm)
{
    using FacetGridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings>;
    bulkFVG.update(facetFVG.gridView(), FacetGridAdapter(gm.getEmbeddings()), true);
}

// main program
int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse command line argument parameters
    Dumux::Parameters::init(argc, argv);

    using BulkGrid = BULKGRIDTYPE;
    using FacetGrid = Dune::FoamGrid<2, 3>;
    using EdgeGrid = Dune::FoamGrid<1, 3>;

    using GridManager = Dumux::FacetCouplingGridManager<BulkGrid, FacetGrid, EdgeGrid>;
    GridManager gridManager;
    gridManager.init();

    // instantiate the grid geometries with caching
    static_assert(!USEBOXINBULK || !USEMPFAINBULK, "Bulk discretization method not uniquely defined");

    using BulkGridView = typename BulkGrid::LeafGridView;
    using TpfaFVGridGeometry = Dumux::CCTpfaFVGridGeometry<BulkGridView, /*caching*/true>;
    using MpfaFVGridGeometry = typename MpfaGridGeometryTraits<BulkGridView>::type;
    using BoxFVGridGeometry = Dumux::BoxFacetCouplingFVGridGeometry<double, BulkGridView, true>;
    using BulkFVGridGeometry = typename std::conditional< USEBOXINBULK,
                                                          BoxFVGridGeometry,
                                                          typename std::conditional<USEMPFAINBULK,
                                                                                    MpfaFVGridGeometry,
                                                                                    TpfaFVGridGeometry>::type >::type;
    BulkFVGridGeometry bulkFvGeometry( gridManager.grid<0>().leafGridView() );

    using FacetGridView = typename FacetGrid::LeafGridView;
    using FacetFVGridGeometry = Dumux::CCTpfaFVGridGeometry<FacetGridView, true>;
    FacetFVGridGeometry facetFvGeometry( gridManager.grid<1>().leafGridView() );
    facetFvGeometry.update();

    using EdgeGridView = typename EdgeGrid::LeafGridView;
    using EdgeFVGridGeometry = Dumux::CCTpfaFVGridGeometry<EdgeGridView, true>;
    EdgeFVGridGeometry edgeFvGeometry( gridManager.grid<2>().leafGridView() );

    // update grid geometries
    edgeFvGeometry.update();
    facetFvGeometry.update();
    updateBulkFvGeometry(bulkFvGeometry, facetFvGeometry, gridManager);

    // instantiate and update mappers for all domain combinations
    Dumux::FacetCouplingMapper<BulkFVGridGeometry, FacetFVGridGeometry> bulkFacetMapper;
    Dumux::FacetCouplingMapper<FacetFVGridGeometry, EdgeFVGridGeometry, /*bulkGridId*/1, /*edgeGridId*/2> facetEdgeMapper;
    Dumux::FacetCouplingThreeDomainMapper<BulkFVGridGeometry, FacetFVGridGeometry, EdgeFVGridGeometry> hierarchyMapper;

    bulkFacetMapper.update(bulkFvGeometry, facetFvGeometry, gridManager.getEmbeddings());
    facetEdgeMapper.update(facetFvGeometry, edgeFvGeometry, gridManager.getEmbeddings());
    hierarchyMapper.update(bulkFvGeometry, facetFvGeometry, edgeFvGeometry, gridManager.getEmbeddings());

    constexpr auto bulkDomainId = Dune::index_constant<0>();
    constexpr auto facetDomainId = Dune::index_constant<1>();
    constexpr auto edgeDomainId = Dune::index_constant<2>();

    // check correctness of bulk-facet map
    for (unsigned int i = 0; i < 2; ++i)
    {
        // check both the map from the bulk facet as well as from the hierarchy mapper
        const auto& bulkFacetMap = i == 0 ? bulkFacetMapper.couplingMap(bulkDomainId, facetDomainId)
                                          : hierarchyMapper.couplingMap(bulkDomainId, facetDomainId);

        const std::size_t expectedNoEntries = USEMPFAINBULK ? 263 : 56;
        if (bulkFacetMap.size() != expectedNoEntries)
            DUNE_THROW(Dune::InvalidStateException, "BulkFacetMap has " << bulkFacetMap.size() << " instead of " << expectedNoEntries << " entries");
        else
            std::cout << "Found " << expectedNoEntries << " entries in bulk-facet map" << std::endl;

        // The rest is too cumbersome to test for mpfa
        if (USEMPFAINBULK)
            continue;

        std::size_t singleCouplings = 0;
        std::size_t doubleCouplings = 0;
        for (const auto& entry : bulkFacetMap)
        {
            const auto cStencilSize = entry.second.couplingStencil.size();

            if (cStencilSize == 1) singleCouplings++;
            else if (cStencilSize == 2) doubleCouplings++;
            else DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 1 or 2");

            const auto bulkElement = bulkFvGeometry.element(entry.first);
            auto fvElementGeometry = localView(bulkFvGeometry);
            fvElementGeometry.bind(bulkElement);

            // check scvf conformity with low dim elements
            for (const auto& elemToScvfs : entry.second.elementToScvfMap)
            {
                const auto lowDimElemIdx = elemToScvfs.first;
                const auto lowDimElement = facetFvGeometry.element(lowDimElemIdx);
                const auto lowDimGeom = lowDimElement.geometry();

                // check correctness of the map
                for (auto scvfIdx : elemToScvfs.second)
                    checkScvfEmbedment(fvElementGeometry.scvf(scvfIdx), lowDimGeom);

                // check correctness of the map using dof indices
                const auto lowDimElemDofIndices = [&] ()
                {
                    std::vector< typename FacetGridView::IndexSet::IndexType > dofIndices;
                    if (FacetFVGridGeometry::discMethod == Dumux::DiscretizationMethod::cctpfa)
                        dofIndices.push_back(lowDimElemIdx);
                    else
                        for (int i = 0; i < lowDimGeom.corners(); ++i)
                            dofIndices.push_back( facetFvGeometry.vertexMapper().subIndex(lowDimElement, i, FacetGridView::dimension) );
                    return dofIndices;
                } ();

                for (auto lowDimDofIdx : lowDimElemDofIndices)
                    for (auto scvfIdx : entry.second.dofToCouplingScvfMap.at(lowDimDofIdx))
                        checkScvfEmbedment(fvElementGeometry.scvf(scvfIdx), lowDimGeom);
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
            const auto lowDimGeom = facetFvGeometry.element(entry.first).geometry();

            // stencil sizes are too cumbersome to test for mpfa
            if (!USEMPFAINBULK)
            {
                const auto cStencilSize = entry.second.couplingStencil.size();
                const std::vector<unsigned int> possibleStencilSizes = USEBOXINBULK ? std::vector<unsigned int>{6, 7, 8}
                                                                                    : std::vector<unsigned int>{2};
                if ( !std::count(possibleStencilSizes.begin(), possibleStencilSizes.end(), cStencilSize) )
                    DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size of " << cStencilSize << " is invalid");
            }

            for (const auto& embedment : entry.second.embedments)
            {
                const auto bulkElement = bulkFvGeometry.element(embedment.first);
                auto fvElementGeometry = localView(bulkFvGeometry);
                fvElementGeometry.bind(bulkElement);

                // check if the scvfs of the embedment coincide with low dim element
                for (auto scvfIdx : embedment.second)
                {
                    const auto& bulkScvf = fvElementGeometry.scvf(scvfIdx);
                    if (fvElementGeometry.scv(bulkScvf.insideScvIdx()).elementIndex() != embedment.first)
                        DUNE_THROW(Dune::InvalidStateException, "Element index in which the scvf is embedded in does not match with bulk element idx.");
                    checkScvfEmbedment(fvElementGeometry.scvf(scvfIdx), lowDimGeom);
                }
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
            const auto bulkElement = facetFvGeometry.element(entry.first);
            auto fvElementGeometry = localView(facetFvGeometry);
            fvElementGeometry.bind(bulkElement);

            const auto cStencilSize = entry.second.couplingStencil.size();
            if (cStencilSize != 1)
                DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 1");

            // we know the edge domain uses tpfa, so coupling stencil = element stencil
            for (unsigned int i = 0; i < cStencilSize; ++i)
            {
                const auto lowDimIdx = entry.second.couplingStencil[i];
                const auto lowDimGeom = edgeFvGeometry.element(lowDimIdx).geometry();

                for (auto scvfIdx : entry.second.dofToCouplingScvfMap.at(lowDimIdx))
                    checkScvfEmbedment(fvElementGeometry.scvf(scvfIdx), lowDimGeom);
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
            const auto lowDimGeom = edgeFvGeometry.element(entry.first).geometry();
            const auto cStencilSize = entry.second.couplingStencil.size();

            if (cStencilSize != 4)
                DUNE_THROW(Dune::InvalidStateException, "Coupling stencil size is " << cStencilSize << " instead of 4");

            for (unsigned int i = 0; i < cStencilSize; ++i)
            {
                const auto& embedment = entry.second.embedments[i];
                const auto facetElement = facetFvGeometry.element(embedment.first);
                auto fvElementGeometry = localView(facetFvGeometry);
                fvElementGeometry.bind(facetElement);

                // check if the scvfs of the embedment coincide with low dim element
                for (auto scvfIdx : embedment.second)
                {
                    const auto& facetScvf = fvElementGeometry.scvf(scvfIdx);
                    if (facetScvf.insideScvIdx() != embedment.first)
                        DUNE_THROW(Dune::InvalidStateException, "Scvf insideScvIdx() does not match with bulk idx from embedment");
                    checkScvfEmbedment(facetScvf, lowDimGeom);
                }
            }
        }
    }

    // everything is ok
    std::cout << "\n... test passed!" << std::endl;
}
