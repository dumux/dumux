// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Pore network of a sphere packing from a regular triangulation of the sphere centres
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_NETWORK_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_NETWORK_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>
#include <dumux/porenetwork/extraction/spherepackingwalls.hh>

namespace Dumux::PoreNetwork::SpherePacking {

/*!
 * \brief Tetrahedra of a triangulation of sphere centres
 *
 * cells[c] holds the vertex indices of tetrahedron c, neighbors[c][j] the tetrahedron sharing
 * the facet opposite to the vertex cells[c][j], or -1 if that facet is on the convex hull.
 * Vertices are spheres, or the walls of a box (see Walls). Hidden spheres have an empty power cell
 * and are not a vertex of any tetrahedron.
 */
struct Triangulation
{
    std::vector<std::array<int, 4>> cells;
    std::vector<std::array<int, 4>> neighbors;
    std::size_t numHiddenSpheres = 0;
};

/*!
 * \brief Pores and throats of a sphere packing
 *
 * For a network of a region, pore labels are -1 in the interior and the label of a boundary face
 * otherwise, throat labels follow from the pore labels as for generated pore networks. For a network
 * bounded by walls, the label of a pore is the sum of 2^k over the walls k it touches, -1 if none, and
 * the label of a throat combines those of its pores.
 */
template<class Scalar>
struct Network
{
    struct Pore
    {
        Point<Scalar> position;
        Scalar volume = 0.0; //!< void volume
        Scalar bulkVolume = 0.0; //!< volume of the tetrahedra including the solid
        Scalar inscribedRadius = 0.0;
        int label = -1;
    };

    struct Throat
    {
        std::array<std::size_t, 2> pores;
        Scalar length = 0.0;
        Scalar fluidArea = 0.0;
        Scalar inscribedRadius = 0.0;
        Scalar hydraulicRadius = 0.0;
        Scalar shapeFactor = 0.0;
        int label = -1;
    };

    //! facet of the triangulation between two different pores, for the pressure forces on its vertices
    struct Facet
    {
        std::array<std::size_t, 2> pores;
        std::array<int, 3> vertices;
        Point<Scalar> normal; //!< unit normal pointing into pores[0]
        std::array<Scalar, 3> forceWeights;
    };

    //! pressure of a pore on a wall it touches: the wall receives the pressure times this vector
    struct WallPressure
    {
        std::size_t pore;
        int vertex;
        Point<Scalar> unitForce;
    };

    std::vector<Pore> pores;
    std::vector<Throat> throats;
    std::vector<Facet> facets;
    std::vector<WallPressure> wallPressures;
    std::vector<long> cellPores; //!< pore of each tetrahedron of the triangulation, -1 if none
    std::size_t numClosedThroats = 0; //!< facets without fluid passage, not part of the network
    std::size_t numMergedTetrahedra = 0; //!< tetrahedra merged into a pore of another tetrahedron
    std::size_t numIsolatedPores = 0; //!< kept tetrahedra without throat, not part of the network
};

/*!
 * \brief Region of a sphere packing turned into a pore network
 *
 * A tetrahedron becomes a pore if its centroid is inside [lower, upper]. It is a boundary pore
 * of the face its vertices reach furthest beyond, labelled with boundaryLabels in the order
 * xMin, xMax, yMin, yMax, zMin, zMax. Tetrahedra whose power centres are closer than
 * mergeDistance form one pore; mergeDistance = 0 merges coinciding power centres only.
 */
template<class Scalar>
struct ExtractionRegion
{
    Point<Scalar> lower;
    Point<Scalar> upper;
    std::array<int, 6> boundaryLabels{1, 2, 3, 4, 5, 6};
    Scalar mergeDistance = 0.0;
};

namespace Detail {

class UnionFind
{
public:
    explicit UnionFind(std::size_t n) : parent_(n)
    { std::iota(parent_.begin(), parent_.end(), 0); }

    std::size_t find(std::size_t i)
    {
        while (parent_[i] != i)
            i = parent_[i] = parent_[parent_[i]];
        return i;
    }

    void unite(std::size_t i, std::size_t j)
    {
        i = find(i); j = find(j);
        if (i != j)
            parent_[std::max(i, j)] = std::min(i, j);
    }

private:
    std::vector<std::size_t> parent_;
};

//! Label of a throat from the labels of its pores, as for generated pore networks
inline int throatLabel(int label0, int label1, const std::array<int, 6>& boundaryLabels)
{
    if (label0 == label1 || label1 == -1)
        return label0;
    if (label0 == -1)
        return label1;
    for (const int label : boundaryLabels)
        if (label0 == label || label1 == label)
            return label;
    return std::max(label0, label1);
}

//! Union of wall masks, -1 for none
inline int combineMasks(int label0, int label1)
{ return std::max(label0, 0) | std::max(label1, 0) ? (std::max(label0, 0) | std::max(label1, 0)) : -1; }

template<class Scalar>
struct CellRecord
{
    bool kept = false;
    typename Network<Scalar>::Pore pore;
    Scalar labelPriority = 0.0;
};

template<class Scalar>
struct FacetRecord
{
    std::array<std::size_t, 2> cells;
    std::array<int, 3> vertices;
    Scalar length, fluidArea, inscribedRadius, hydraulicRadius, shapeFactor;
    Point<Scalar> normal;
    std::array<Scalar, 3> forceWeights;
};

template<class Scalar>
struct WallRecord
{
    std::size_t cell;
    int vertex;
    Point<Scalar> unitForce;
};

/*!
 * \brief Network from the pores of single tetrahedra and the facets between them
 *
 * Tetrahedra joined by a facet shorter than mergeDistance form one pore. Parallel facets between two
 * pores are combined into one throat with the summed fluid area, the area-weighted mean length and the
 * hydraulic radius that conserves the sum of A R_h^2/L. Pores without throat are left out.
 */
template<class Scalar>
Network<Scalar> assembleNetwork(const std::vector<CellRecord<Scalar>>& cells,
                                const std::vector<FacetRecord<Scalar>>& facets,
                                const std::vector<WallRecord<Scalar>>& wallRecords,
                                Scalar mergeDistance, bool labelsAreMasks,
                                const std::array<int, 6>& boundaryLabels)
{
    using Pore = typename Network<Scalar>::Pore;
    const std::size_t numCells = cells.size();
    Network<Scalar> network;

    UnionFind groups(numCells);
    for (const auto& f : facets)
        if (f.length <= mergeDistance)
            groups.unite(f.cells[0], f.cells[1]);

    // pores of the groups of merged tetrahedra
    std::vector<std::size_t> groupOf(numCells);
    std::map<std::size_t, std::size_t> groupIndex;
    std::vector<std::size_t> groupSize;
    std::vector<Pore> groupPores;
    std::vector<Scalar> groupPriority;
    for (std::size_t c = 0; c < numCells; ++c)
    {
        if (!cells[c].kept)
            continue;
        const auto& cellPore = cells[c].pore;
        auto [it, inserted] = groupIndex.try_emplace(groups.find(c), groupPores.size());
        if (inserted)
        {
            groupPores.push_back(cellPore);
            groupPores.back().position = 0.0;
            groupPores.back().volume = 0.0;
            groupPores.back().bulkVolume = 0.0;
            groupSize.push_back(0);
            groupPriority.push_back(cells[c].labelPriority);
        }
        const auto g = it->second;
        groupOf[c] = g;
        auto& pore = groupPores[g];
        pore.position += cellPore.position;
        pore.volume += cellPore.volume;
        pore.bulkVolume += cellPore.bulkVolume;
        pore.inscribedRadius = std::max(pore.inscribedRadius, cellPore.inscribedRadius);
        if (labelsAreMasks)
            pore.label = combineMasks(pore.label, cellPore.label);
        else if (cells[c].labelPriority > groupPriority[g])
        {
            groupPriority[g] = cells[c].labelPriority;
            pore.label = cellPore.label;
        }
        ++groupSize[g];
    }
    for (std::size_t g = 0; g < groupPores.size(); ++g)
    {
        groupPores[g].position /= static_cast<Scalar>(groupSize[g]);
        network.numMergedTetrahedra += groupSize[g] - 1;
    }

    // throats between groups, parallel facets combined
    struct Combined
    {
        std::size_t count = 0;
        const FacetRecord<Scalar>* single = nullptr;
        Scalar area = 0.0, areaLength = 0.0, conductanceFactor = 0.0, areaShapeFactor = 0.0, inscribedRadius = 0.0;
    };
    std::map<std::pair<std::size_t, std::size_t>, Combined> combined;
    for (const auto& f : facets)
    {
        const auto g0 = groupOf[f.cells[0]], g1 = groupOf[f.cells[1]];
        if (g0 == g1)
            continue;
        if (f.inscribedRadius <= 0.0 || f.hydraulicRadius <= 0.0)
        {
            ++network.numClosedThroats;
            continue;
        }
        auto& entry = combined[std::minmax(g0, g1)];
        ++entry.count;
        entry.single = &f;
        entry.area += f.fluidArea;
        entry.areaLength += f.fluidArea*f.length;
        entry.conductanceFactor += f.fluidArea*f.hydraulicRadius*f.hydraulicRadius/f.length;
        entry.areaShapeFactor += f.fluidArea*f.shapeFactor;
        entry.inscribedRadius = std::max(entry.inscribedRadius, f.inscribedRadius);
    }

    // only pores with at least one throat are part of the network
    std::vector<long> poreIndex(groupPores.size(), -1);
    for (const auto& [pair, entry] : combined)
        poreIndex[pair.first] = poreIndex[pair.second] = 0;
    for (std::size_t g = 0; g < groupPores.size(); ++g)
    {
        if (poreIndex[g] < 0)
        {
            ++network.numIsolatedPores;
            continue;
        }
        poreIndex[g] = network.pores.size();
        network.pores.push_back(groupPores[g]);
    }

    for (const auto& [pair, entry] : combined)
    {
        typename Network<Scalar>::Throat throat;
        throat.pores = {static_cast<std::size_t>(poreIndex[pair.first]), static_cast<std::size_t>(poreIndex[pair.second])};
        if (entry.count == 1)
        {
            throat.length = entry.single->length;
            throat.fluidArea = entry.single->fluidArea;
            throat.inscribedRadius = entry.single->inscribedRadius;
            throat.hydraulicRadius = entry.single->hydraulicRadius;
            throat.shapeFactor = entry.single->shapeFactor;
        }
        else
        {
            using std::sqrt;
            throat.length = entry.areaLength/entry.area;
            throat.fluidArea = entry.area;
            throat.inscribedRadius = entry.inscribedRadius;
            throat.hydraulicRadius = sqrt(entry.conductanceFactor*throat.length/entry.area);
            throat.shapeFactor = entry.areaShapeFactor/entry.area;
        }
        const int label0 = network.pores[throat.pores[0]].label, label1 = network.pores[throat.pores[1]].label;
        throat.label = labelsAreMasks ? combineMasks(label0, label1) : throatLabel(label0, label1, boundaryLabels);
        network.throats.push_back(throat);
    }

    for (const auto& f : facets)
    {
        const auto g0 = groupOf[f.cells[0]], g1 = groupOf[f.cells[1]];
        if (g0 == g1 || poreIndex[g0] < 0 || poreIndex[g1] < 0)
            continue;
        network.facets.push_back({{static_cast<std::size_t>(poreIndex[g0]), static_cast<std::size_t>(poreIndex[g1])},
                                  f.vertices, f.normal, f.forceWeights});
    }

    network.cellPores.assign(numCells, -1);
    for (std::size_t c = 0; c < numCells; ++c)
        if (cells[c].kept && poreIndex[groupOf[c]] >= 0)
            network.cellPores[c] = poreIndex[groupOf[c]];

    for (const auto& w : wallRecords)
        if (cells[w.cell].kept && poreIndex[groupOf[w.cell]] >= 0)
            network.wallPressures.push_back({static_cast<std::size_t>(poreIndex[groupOf[w.cell]]), w.vertex, w.unitForce});

    return network;
}

} // end namespace Detail

/*!
 * \brief Pore network of the tetrahedra of a regular triangulation of the spheres in a region
 *
 * Pores: power centre, void and bulk volume and inscribed sphere of the tetrahedron; for a flat
 * tetrahedron without a sphere touching all four spheres, the largest sphere around the power centre.
 * Throats: facets shared by two pores, with the distance of the power centres as length, the facet
 * fluid area, the inscribed circle, the hydraulic radius of the region between facet and power
 * centres, and the shape factor A/P^2 with the wetted perimeter P of the facet.
 * Parallel throats between merged pores are combined into one throat with the summed fluid area,
 * the area-weighted mean length and the hydraulic radius that conserves the sum of A R_h^2/L.
 */
template<class Scalar>
Network<Scalar> extractNetwork(const Triangulation& triangulation,
                               const std::vector<Point<Scalar>>& centers,
                               const std::vector<Scalar>& radii,
                               const ExtractionRegion<Scalar>& region)
{
    const std::size_t numCells = triangulation.cells.size();
    std::vector<Detail::CellRecord<Scalar>> cells(numCells);
    for (std::size_t c = 0; c < numCells; ++c)
    {
        std::array<Point<Scalar>, 4> x;
        std::array<Scalar, 4> r;
        for (int j = 0; j < 4; ++j)
        {
            x[j] = centers[triangulation.cells[c][j]];
            r[j] = radii[triangulation.cells[c][j]];
        }
        const auto centroid = (x[0] + x[1] + x[2] + x[3])/4.0;
        bool inside = true;
        for (int axis = 0; axis < 3; ++axis)
            inside = inside && centroid[axis] >= region.lower[axis] && centroid[axis] <= region.upper[axis];
        if (!inside)
            continue;

        auto& cell = cells[c];
        cell.kept = true;
        cell.pore.position = powerCenter(x, r);
        cell.pore.bulkVolume = tetrahedronVolume(x);
        cell.pore.volume = cell.pore.bulkVolume - tetrahedronSolidVolume(x, r);
        cell.pore.inscribedRadius = inscribedSphere(x, r).radius;
        if (cell.pore.inscribedRadius <= 0.0)
            cell.pore.inscribedRadius = clearance(cell.pore.position, x, r);
        for (int face = 0; face < 6; ++face)
        {
            const int axis = face/2;
            for (const auto& v : x)
            {
                const Scalar depth = face % 2 == 0 ? region.lower[axis] - v[axis] : v[axis] - region.upper[axis];
                if (depth > cell.labelPriority)
                {
                    cell.labelPriority = depth;
                    cell.pore.label = region.boundaryLabels[face];
                }
            }
        }
    }

    std::vector<Detail::FacetRecord<Scalar>> facets;
    for (std::size_t c = 0; c < numCells; ++c)
    {
        if (!cells[c].kept)
            continue;
        for (int j = 0; j < 4; ++j)
        {
            const int n = triangulation.neighbors[c][j];
            if (n < 0 || !cells[n].kept || static_cast<std::size_t>(n) < c)
                continue;

            Detail::FacetRecord<Scalar> facet;
            facet.cells = {c, static_cast<std::size_t>(n)};
            std::array<Point<Scalar>, 3> x;
            std::array<Scalar, 3> r;
            for (int k = 0, m = 0; k < 4; ++k)
            {
                if (k == j)
                    continue;
                facet.vertices[m] = triangulation.cells[c][k];
                x[m] = centers[facet.vertices[m]];
                r[m] = radii[facet.vertices[m]];
                ++m;
            }

            const auto& p1 = cells[c].pore.position;
            const auto& p2 = cells[n].pore.position;
            facet.length = (p1 - p2).two_norm();
            facet.fluidArea = facetFluidArea(x, r);
            facet.inscribedRadius = inscribedCircle(x, r).radius;
            facet.hydraulicRadius = throatRegion(x, r, p1, p2).hydraulicRadius();
            Scalar wettedPerimeter = 0.0;
            for (int k = 0; k < 3; ++k)
                wettedPerimeter += r[k]*planeAngle(x[k], x[(k+1)%3], x[(k+2)%3]);
            facet.shapeFactor = facet.fluidArea/(wettedPerimeter*wettedPerimeter);
            const auto force = facetForce(x, r, p1, p2);
            facet.normal = force.normal;
            facet.forceWeights = force.weights;
            facets.push_back(facet);
        }
    }

    return Detail::assembleNetwork(cells, facets, std::vector<Detail::WallRecord<Scalar>>{},
                                   region.mergeDistance, false, region.boundaryLabels);
}

/*!
 * \brief Pore network of a sphere packing in a box, triangulated together with the walls of the box
 *        (see regularTriangulation with walls)
 *
 * All tetrahedra are pores, those touching walls with the geometry of spherepackingwalls.hh; facets of
 * three walls lie on the convex hull. Pore labels are wall masks (see Network). The forces on the walls
 * are the pressures of the pores on the projections of their facets opposite to the walls.
 */
template<class Scalar>
Network<Scalar> extractNetwork(const Triangulation& triangulation,
                               const std::vector<Point<Scalar>>& centers,
                               const std::vector<Scalar>& radii,
                               const Walls<Scalar>& walls,
                               Scalar mergeDistance = 0.0)
{
    const std::size_t numCells = triangulation.cells.size();
    const auto vertex = [&](int v) {
        return walls.isWall(v) ? Vertex<Scalar>{walls.center(walls.wall(v)), walls.radius(), walls.wall(v)}
                               : Vertex<Scalar>{centers[v], radii[v], -1};
    };
    const auto cellVertices = [&](std::size_t c) {
        std::array<Vertex<Scalar>, 4> v;
        for (int j = 0; j < 4; ++j)
            v[j] = vertex(triangulation.cells[c][j]);
        return v;
    };

    // power centre with a sphere as first vertex, which keeps the far wall centres out of the shift
    const auto cellPowerCenter = [&](const std::array<Vertex<Scalar>, 4>& v) {
        int first = 0;
        while (v[first].wall >= 0)
            ++first;
        std::array<Point<Scalar>, 4> x;
        std::array<Scalar, 4> r;
        for (int j = 0; j < 4; ++j)
        {
            x[j] = v[(first + j)%4].x;
            r[j] = v[(first + j)%4].r;
        }
        return powerCenter(x, r);
    };

    std::vector<Detail::CellRecord<Scalar>> cells(numCells);
    std::vector<Point<Scalar>> powerCenters(numCells);
    for (std::size_t c = 0; c < numCells; ++c)
    {
        const auto v = cellVertices(c);
        const int numWalls = std::count_if(v.begin(), v.end(), [](const auto& vertex) { return vertex.wall >= 0; });
        if (numWalls == 4)
            continue;

        powerCenters[c] = cellPowerCenter(v);
        auto& cell = cells[c];
        cell.kept = true;
        cell.pore.position = powerCenters[c];
        cell.pore.bulkVolume = cellBulkVolume(v, walls);
        cell.pore.volume = cell.pore.bulkVolume - cellSolidVolume(v);
        if (numWalls == 0)
            cell.pore.inscribedRadius = inscribedSphere(std::array<Point<Scalar>, 4>{v[0].x, v[1].x, v[2].x, v[3].x},
                                                        std::array<Scalar, 4>{v[0].r, v[1].r, v[2].r, v[3].r}).radius;
        else
            cell.pore.inscribedRadius = cellInscribedSphere(v, walls).radius;
        if (cell.pore.inscribedRadius <= 0.0)
        {
            std::vector<Point<Scalar>> x;
            std::vector<Scalar> r;
            for (const auto& vertex : v)
                if (vertex.wall < 0)
                {
                    x.push_back(vertex.x);
                    r.push_back(vertex.r);
                }
            cell.pore.inscribedRadius = clearance(powerCenters[c], x, r);
        }
        int mask = 0;
        for (const auto& vertex : v)
            if (vertex.wall >= 0)
                mask |= 1 << vertex.wall;
        cell.pore.label = mask > 0 ? mask : -1;
    }

    const auto facetVertices = [&](std::size_t c, int j) {
        std::array<Vertex<Scalar>, 3> v;
        std::array<int, 3> ids;
        for (int k = 0, m = 0; k < 4; ++k)
            if (k != j)
            {
                ids[m] = triangulation.cells[c][k];
                v[m] = vertex(ids[m]);
                ++m;
            }
        return std::make_pair(v, ids);
    };

    std::vector<Detail::FacetRecord<Scalar>> facets;
    std::vector<Detail::WallRecord<Scalar>> wallRecords;
    for (std::size_t c = 0; c < numCells; ++c)
    {
        if (!cells[c].kept)
            continue;
        for (int j = 0; j < 4; ++j)
        {
            const int n = triangulation.neighbors[c][j];
            if (n < 0 || !cells[n].kept)
                continue;
            const auto [v, ids] = facetVertices(c, j);
            const int numWalls = std::count_if(v.begin(), v.end(), [](const auto& vertex) { return vertex.wall >= 0; });
            if (numWalls == 3)
                continue;
            const auto& p1 = powerCenters[c];
            const auto& p2 = powerCenters[n];
            const auto geometry = facetGeometry(v, p1, p2, walls);

            // pressure on the wall opposite to the facet, over the projection of the facet onto it
            const int opposite = triangulation.cells[c][j];
            if (walls.isWall(opposite))
            {
                using std::abs;
                const int w = walls.wall(opposite);
                wallRecords.push_back({c, opposite, -abs(geometry.areaVector[Walls<Scalar>::axis(w)])*walls.inwardNormal(w)});
            }

            if (static_cast<std::size_t>(n) < c)
                continue;

            Detail::FacetRecord<Scalar> facet;
            facet.cells = {c, static_cast<std::size_t>(n)};
            facet.vertices = ids;
            facet.length = (p1 - p2).two_norm();
            facet.fluidArea = geometry.fluidArea;
            facet.inscribedRadius = geometry.inscribedRadius;
            facet.hydraulicRadius = geometry.hydraulicRadius;
            Scalar wettedPerimeter = 0.0, solidSurface = 0.0;
            for (int k = 0; k < 3; ++k)
            {
                if (v[k].wall < 0)
                    wettedPerimeter += 2.0*geometry.crossSection[k]/v[k].r;
                solidSurface += geometry.solidSurface[k];
            }
            facet.shapeFactor = facet.fluidArea/(wettedPerimeter*wettedPerimeter);
            facet.normal = geometry.areaVector/geometry.areaVector.two_norm();
            for (int k = 0; k < 3; ++k)
                facet.forceWeights[k] = geometry.crossSection[k] + facet.fluidArea*geometry.solidSurface[k]/solidSurface;
            facets.push_back(facet);
        }
    }

    return Detail::assembleNetwork(cells, facets, wallRecords, mergeDistance, true, std::array<int, 6>{});
}

/*!
 * \brief Bulk volumes of the pores for new positions of the spheres and walls, with the tetrahedra of the
 *        triangulation the network was extracted from
 *
 * The rate of change of these volumes drives the flow when the packing deforms (Catalano et al. 2014).
 */
template<class Scalar>
std::vector<Scalar> poreBulkVolumes(const Network<Scalar>& network, const Triangulation& triangulation,
                                    const std::vector<Point<Scalar>>& centers, const std::vector<Scalar>& radii,
                                    const Walls<Scalar>* walls = nullptr)
{
    std::vector<Scalar> volumes(network.pores.size(), 0.0);
    for (std::size_t c = 0; c < triangulation.cells.size(); ++c)
    {
        const long pore = network.cellPores[c];
        if (pore < 0)
            continue;
        std::array<Vertex<Scalar>, 4> v;
        for (int j = 0; j < 4; ++j)
        {
            const int id = triangulation.cells[c][j];
            v[j] = walls && walls->isWall(id) ? Vertex<Scalar>{walls->center(walls->wall(id)), walls->radius(), walls->wall(id)}
                                              : Vertex<Scalar>{centers[id], radii[id], -1};
        }
        volumes[pore] += walls ? cellBulkVolume(v, *walls)
                               : tetrahedronVolume(std::array<Point<Scalar>, 4>{v[0].x, v[1].x, v[2].x, v[3].x});
    }
    return volumes;
}

/*!
 * \brief Pressure forces of the pores on the spheres and walls (Chareyre et al. 2012, without viscous shear)
 *
 * Every facet between two pores adds -(p_0 - p_1) n w_k to each of its vertices k (see facetForce), and
 * every pore touching a wall pushes it with its pressure on the projection of the facet opposite to the
 * wall. Forces are returned for numVertices vertices: the spheres, followed by the walls if present.
 * Spheres next to the boundary of a network region only receive the forces of the facets inside it.
 */
template<class Scalar, class PressureVector>
std::vector<Point<Scalar>> fluidForces(const Network<Scalar>& network, const PressureVector& pressure,
                                       std::size_t numVertices)
{
    std::vector<Point<Scalar>> forces(numVertices, Point<Scalar>(0.0));
    for (const auto& f : network.facets)
    {
        const Scalar dp = pressure[f.pores[0]] - pressure[f.pores[1]];
        for (int k = 0; k < 3; ++k)
            forces[f.vertices[k]].axpy(-dp*f.forceWeights[k], f.normal);
    }
    for (const auto& w : network.wallPressures)
        forces[w.vertex].axpy(pressure[w.pore], w.unitForce);
    return forces;
}

/*!
 * \brief Write the network as DGF file readable by PoreNetwork::GridManager
 *
 * Pore parameters: PoreInscribedRadius PoreVolume PoreLabel;
 * throat parameters: ThroatInscribedRadius ThroatLength ThroatCrossSectionalArea ThroatShapeFactor
 * ThroatHydraulicRadius ThroatLabel. The pore geometry is not written, set Grid.PoreGeometry.
 */
template<class Scalar>
void writeDgf(const std::string& fileName, const Network<Scalar>& network)
{
    std::ofstream file(fileName);
    if (!file)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

    file << std::setprecision(17);
    file << "DGF\n"
         << "% Vertex parameters: PoreInscribedRadius PoreVolume PoreLabel\n"
         << "% Element parameters: ThroatInscribedRadius ThroatLength ThroatCrossSectionalArea ThroatShapeFactor ThroatHydraulicRadius ThroatLabel\n"
         << "Vertex % pores of a sphere packing\n"
         << "parameters 3\n";
    for (const auto& pore : network.pores)
        file << pore.position[0] << " " << pore.position[1] << " " << pore.position[2] << " "
             << pore.inscribedRadius << " " << pore.volume << " " << pore.label << "\n";
    file << "#\n"
         << "SIMPLEX % throats of a sphere packing\n"
         << "parameters 6\n";
    for (const auto& throat : network.throats)
        file << throat.pores[0] << " " << throat.pores[1] << " "
             << throat.inscribedRadius << " " << throat.length << " " << throat.fluidArea << " "
             << throat.shapeFactor << " " << throat.hydraulicRadius << " " << throat.label << "\n";
    file << "#\n";
}

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
