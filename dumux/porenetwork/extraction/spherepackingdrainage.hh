// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Quasi-static drainage of the pore network of a sphere packing
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_DRAINAGE_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_DRAINAGE_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/porenetwork/2p/static/staticdrainge.hh>
#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>
#include <dumux/porenetwork/extraction/spherepackingnetwork.hh>
#include <dumux/porenetwork/extraction/spherepackingwalls.hh>

namespace Dumux::PoreNetwork::SpherePacking {

/*!
 * \brief Drainage entry radius of every throat, the entry capillary pressure is surfaceTension over it
 *
 * Facets between three spheres take the Mayer-Stowe-Princen radius (mspThroatRadius), facets with
 * walls the radius of their inscribed circle, as YADE's TwoPhaseFlowEngine (entryPressureMethod 1).
 * A throat combining several facets takes the largest radius of its facets.
 */
template<class Scalar>
std::vector<Scalar> throatEntryRadii(const Network<Scalar>& network,
                                     const std::vector<Point<Scalar>>& centers,
                                     const std::vector<Scalar>& radii,
                                     const Walls<Scalar>& walls)
{
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> throatIndex;
    for (std::size_t t = 0; t < network.throats.size(); ++t)
        throatIndex[std::minmax(network.throats[t].pores[0], network.throats[t].pores[1])] = t;

    std::vector<Scalar> entryRadius(network.throats.size(), 0.0);
    for (const auto& facet : network.facets)
    {
        std::array<Vertex<Scalar>, 3> v;
        int numWalls = 0;
        for (int k = 0; k < 3; ++k)
        {
            const int id = facet.vertices[k];
            v[k] = walls.isWall(id) ? Vertex<Scalar>{walls.center(walls.wall(id)), walls.radius(), walls.wall(id)}
                                    : Vertex<Scalar>{centers[id], radii[id], -1};
            numWalls += walls.isWall(id);
        }
        const Scalar radius = numWalls == 0
            ? mspThroatRadius(std::array<Point<Scalar>, 3>{v[0].x, v[1].x, v[2].x}, std::array<Scalar, 3>{v[0].r, v[1].r, v[2].r})
            : facetGeometry(v, network.pores[facet.pores[0]].position, network.pores[facet.pores[1]].position, walls).inscribedRadius;
        const auto it = throatIndex.find(std::minmax(facet.pores[0], facet.pores[1]));
        if (it != throatIndex.end())
            entryRadius[it->second] = std::max(entryRadius[it->second], radius);
    }
    return entryRadius;
}

/*!
 * \brief Quasi-static drainage of the pore network of a sphere packing at a given capillary pressure,
 *        computed with TwoPStaticDrainage on a one-dimensional grid of the network
 *
 * Inlet pores belong to the non-wetting reservoir and are drained at any capillary pressure, outlet
 * pores belong to the wetting reservoir and are never drained. Another pore drains when the capillary
 * pressure reaches the entry pressure of a throat to a drained pore. The wetting phase is not trapped.
 */
template<class Scalar>
class NetworkDrainage
{
    using Grid = Dune::FoamGrid<1, 3>;

    struct GridGeometry
    {
        using GridView = typename Grid::LeafGridView;
        explicit GridGeometry(const GridView& gv) : gridView_(gv) {}
        const GridView& gridView() const { return gridView_; }
    private:
        GridView gridView_;
    };

public:
    enum ThroatLabel { interior = 0, inlet = 1, outlet = 2 };

    NetworkDrainage(const Network<Scalar>& network,
                    const std::vector<Scalar>& throatEntryPressures,
                    const std::vector<bool>& isInletPore,
                    const std::vector<bool>& isOutletPore)
    : isInletPore_(isInletPore)
    , isOutletPore_(isOutletPore)
    {
        if (throatEntryPressures.size() != network.throats.size()
            || isInletPore.size() != network.pores.size() || isOutletPore.size() != network.pores.size())
            DUNE_THROW(Dune::InvalidStateException, "Sizes of the drainage data do not match the network");

        Dune::GridFactory<Grid> factory;
        for (const auto& pore : network.pores)
            factory.insertVertex(pore.position);
        for (const auto& throat : network.throats)
            factory.insertElement(Dune::GeometryTypes::line,
                                  {static_cast<unsigned int>(throat.pores[0]), static_cast<unsigned int>(throat.pores[1])});
        grid_ = factory.createGrid();
        gridGeometry_ = std::make_unique<GridGeometry>(grid_->leafGridView());

        const auto& gridView = gridGeometry_->gridView();
        elementThroat_.resize(gridView.size(0));
        entryPressure_.resize(gridView.size(0));
        label_.resize(gridView.size(0));
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(element);
            const auto t = factory.insertionIndex(element);
            const auto& pores = network.throats[t].pores;
            elementThroat_[eIdx] = t;
            entryPressure_[eIdx] = throatEntryPressures[t];
            label_[eIdx] = isInletPore[pores[0]] || isInletPore[pores[1]] ? inlet
                         : isOutletPore[pores[0]] || isOutletPore[pores[1]] ? outlet : interior;
        }
        throatPores_.reserve(network.throats.size());
        for (const auto& throat : network.throats)
            throatPores_.push_back(throat.pores);
    }

    //! drained state of every pore at the capillary pressure pc
    std::vector<bool> drainedPores(Scalar pc) const
    {
        TwoPStaticDrainage<GridGeometry, Scalar> drainage(*gridGeometry_, entryPressure_, label_, inlet, outlet);
        std::vector<bool> invaded(entryPressure_.size(), false);
        drainage.updateInvasionState(invaded, pc);

        std::vector<bool> drained = isInletPore_;
        for (std::size_t e = 0; e < invaded.size(); ++e)
            if (invaded[e])
                for (const auto pore : throatPores_[elementThroat_[e]])
                    drained[pore] = true;
        for (std::size_t p = 0; p < drained.size(); ++p)
            if (isOutletPore_[p])
                drained[p] = false;
        return drained;
    }

private:
    std::vector<bool> isInletPore_;
    std::vector<bool> isOutletPore_;
    std::unique_ptr<Grid> grid_;
    std::unique_ptr<GridGeometry> gridGeometry_;
    std::vector<std::size_t> elementThroat_;
    std::vector<Scalar> entryPressure_;
    std::vector<int> label_;
    std::vector<std::array<std::size_t, 2>> throatPores_;
};

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
