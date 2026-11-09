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
#include <map>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

#include <dumux/io/grid/porenetwork/networkdata.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>

namespace Dumux::PoreNetwork::SpherePacking {

/*!
 * \brief Tetrahedra of a triangulation of sphere centres
 *
 * cells[c] holds the vertex indices of tetrahedron c, neighbors[c][j] the tetrahedron sharing
 * the facet opposite to the vertex cells[c][j], or -1 if that facet is on the convex hull.
 * Vertices are the indices of the spheres. Hidden spheres have an empty power cell and are not a
 * vertex of any tetrahedron.
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
 * Pore labels are -1 in the interior and the label of a boundary face otherwise, throat labels follow
 * from the pore labels as for generated pore networks (see assembleNetwork).
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
        std::optional<Scalar> transmissibility; //!< at unit viscosity if given, see throatTransmissibility
    };

    //! facet of the triangulation between two different pores
    struct Facet
    {
        std::array<std::size_t, 2> pores;
        std::array<std::size_t, 2> cells; //!< the tetrahedra on the sides of pores[0] and pores[1]
        std::array<int, 3> vertices;
    };

    std::vector<Pore> pores;
    std::vector<Throat> throats;
    std::vector<Facet> facets;
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

//! Pore of one tetrahedron of a triangulation, input of assembleNetwork
template<class Scalar>
struct TetrahedronPore
{
    bool kept = false; //!< whether the tetrahedron is part of the network
    typename Network<Scalar>::Pore pore;
    Scalar labelPriority = 0.0; //!< among merged tetrahedra, the label of the highest priority wins
};

//! Throat of one facet between two tetrahedra of a triangulation, input of assembleNetwork
template<class Scalar>
struct FacetThroat
{
    std::array<std::size_t, 2> cells;
    std::array<int, 3> vertices;
    Scalar length, fluidArea, inscribedRadius, hydraulicRadius, shapeFactor;
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

} // end namespace Detail

/*!
 * \brief Network from the pores of single tetrahedra and the facets between them
 *
 * Tetrahedra joined by a facet shorter than mergeDistance form one pore at the mean of their positions,
 * with the summed volumes, the largest inscribed radius and the label of the tetrahedron with the highest
 * label priority. Facets without fluid passage (zero inscribed or hydraulic radius) are closed. Parallel
 * facets between two pores are combined into one throat with the summed fluid area, the area-weighted mean
 * length and shape factor, and the hydraulic radius that conserves the sum of A R_h^2/L. Pores without
 * throat are left out. Throat labels follow from the pore labels as for generated pore networks, with the
 * labels of the boundary faces in boundaryLabels.
 */
template<class Scalar>
Network<Scalar> assembleNetwork(const std::vector<TetrahedronPore<Scalar>>& cells,
                                const std::vector<FacetThroat<Scalar>>& facets,
                                Scalar mergeDistance,
                                const std::array<int, 6>& boundaryLabels)
{
    using Pore = typename Network<Scalar>::Pore;
    const std::size_t numCells = cells.size();
    Network<Scalar> network;

    Detail::UnionFind groups(numCells);
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
        if (cells[c].labelPriority > groupPriority[g])
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
        const FacetThroat<Scalar>* single = nullptr;
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
        throat.label = Detail::throatLabel(label0, label1, boundaryLabels);
        network.throats.push_back(throat);
    }

    for (const auto& f : facets)
    {
        const auto g0 = groupOf[f.cells[0]], g1 = groupOf[f.cells[1]];
        if (g0 == g1 || poreIndex[g0] < 0 || poreIndex[g1] < 0)
            continue;
        network.facets.push_back({{static_cast<std::size_t>(poreIndex[g0]), static_cast<std::size_t>(poreIndex[g1])},
                                  f.cells, f.vertices});
    }

    network.cellPores.assign(numCells, -1);
    for (std::size_t c = 0; c < numCells; ++c)
        if (cells[c].kept && poreIndex[groupOf[c]] >= 0)
            network.cellPores[c] = poreIndex[groupOf[c]];

    return network;
}

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
    std::vector<TetrahedronPore<Scalar>> cells(numCells);
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

    std::vector<FacetThroat<Scalar>> facets;
    for (std::size_t c = 0; c < numCells; ++c)
    {
        if (!cells[c].kept)
            continue;
        for (int j = 0; j < 4; ++j)
        {
            const int n = triangulation.neighbors[c][j];
            if (n < 0 || !cells[n].kept || static_cast<std::size_t>(n) < c)
                continue;

            FacetThroat<Scalar> facet;
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
            facets.push_back(facet);
        }
    }

    return assembleNetwork(cells, facets, region.mergeDistance, region.boundaryLabels);
}

/*!
 * \brief Single-phase transmissibility of a throat at unit viscosity: the given one, otherwise
 *        TransmissibilityChareyre of the throat geometry
 */
template<class Throat>
auto throatTransmissibility(const Throat& throat)
{
    using Scalar = decltype(throat.length);
    return throat.transmissibility.value_or(
        TransmissibilityChareyre<Scalar>::singlePhaseTransmissibility(throat.fluidArea, throat.hydraulicRadius, throat.length));
}

/*!
 * \brief The pores and throats of the network with their parameters, for PoreNetwork::GridManager::init
 *
 * Pore parameters: PoreInscribedRadius PoreVolume PoreLabel;
 * throat parameters: ThroatInscribedRadius ThroatLength ThroatCrossSectionalArea ThroatShapeFactor
 * ThroatHydraulicRadius ThroatLabel ThroatTransmissibility (see throatTransmissibility).
 * The pore geometry is not included, set Grid.PoreGeometry.
 */
template<class Scalar>
PoreNetwork::NetworkData<3, Scalar> networkData(const Network<Scalar>& network)
{
    PoreNetwork::NetworkData<3, Scalar> data;
    data.vertexParameterNames = {"PoreInscribedRadius", "PoreVolume", "PoreLabel"};
    data.elementParameterNames = {"ThroatInscribedRadius", "ThroatLength", "ThroatCrossSectionalArea",
                                  "ThroatShapeFactor", "ThroatHydraulicRadius", "ThroatLabel", "ThroatTransmissibility"};
    for (const auto& pore : network.pores)
    {
        data.porePositions.push_back(pore.position);
        data.vertexParameters.push_back({pore.inscribedRadius, pore.volume, static_cast<Scalar>(pore.label)});
    }
    for (const auto& throat : network.throats)
    {
        data.throatPores.push_back(throat.pores);
        data.elementParameters.push_back({throat.inscribedRadius, throat.length, throat.fluidArea, throat.shapeFactor,
                                          throat.hydraulicRadius, static_cast<Scalar>(throat.label),
                                          throatTransmissibility(throat)});
    }
    return data;
}

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
