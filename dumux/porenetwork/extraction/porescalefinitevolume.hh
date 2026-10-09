// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Fluid coupling of a deforming sphere packing in a box (pore-scale finite volume method; requires
 *        CGAL and UMFPack)
 */
#ifndef DUMUX_PNM_EXTRACTION_PORE_SCALE_FINITE_VOLUME_HH
#define DUMUX_PNM_EXTRACTION_PORE_SCALE_FINITE_VOLUME_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/porenetwork/extraction/porescaleflow.hh>
#include <dumux/porenetwork/extraction/spherepackingnetwork.hh>
#include <dumux/porenetwork/extraction/spherepackingtriangulation.hh>

namespace Dumux::PoreNetwork::SpherePacking {

//! Index of the nearest reference point for each point, searched in a grid of buckets around the reference points
template<class Scalar>
std::vector<std::size_t> nearestPoints(const std::vector<Point<Scalar>>& reference, const std::vector<Point<Scalar>>& points)
{
    using std::floor; using std::cbrt; using std::max;
    Point<Scalar> lower(std::numeric_limits<Scalar>::max()), upper(std::numeric_limits<Scalar>::lowest());
    for (const auto& x : reference)
        for (int c = 0; c < 3; ++c)
        {
            lower[c] = std::min(lower[c], x[c]);
            upper[c] = max(upper[c], x[c]);
        }
    const auto extent = upper - lower;
    const Scalar h = max(cbrt(max(extent[0]*extent[1]*extent[2], Scalar(1e-300))/reference.size()), Scalar(1e-300));
    std::array<long, 3> n;
    for (int c = 0; c < 3; ++c)
        n[c] = static_cast<long>(extent[c]/h) + 1;
    const auto bucketOf = [&](const Point<Scalar>& x, int c) {
        return std::clamp(static_cast<long>(floor((x[c] - lower[c])/h)), 0L, n[c] - 1);
    };
    std::vector<std::vector<std::size_t>> buckets(n[0]*n[1]*n[2]);
    for (std::size_t i = 0; i < reference.size(); ++i)
        buckets[(bucketOf(reference[i], 2)*n[1] + bucketOf(reference[i], 1))*n[0] + bucketOf(reference[i], 0)].push_back(i);

    std::vector<std::size_t> nearest(points.size());
    for (std::size_t p = 0; p < points.size(); ++p)
    {
        const auto& x = points[p];
        const std::array<long, 3> b{bucketOf(x, 0), bucketOf(x, 1), bucketOf(x, 2)};
        Scalar best = std::numeric_limits<Scalar>::max();
        for (long ring = 0; ; ++ring)
        {
            for (long k = b[2] - ring; k <= b[2] + ring; ++k)
                for (long j = b[1] - ring; j <= b[1] + ring; ++j)
                    for (long i = b[0] - ring; i <= b[0] + ring; ++i)
                    {
                        if (std::max({std::abs(i - b[0]), std::abs(j - b[1]), std::abs(k - b[2])}) != ring
                            || i < 0 || j < 0 || k < 0 || i >= n[0] || j >= n[1] || k >= n[2])
                            continue;
                        for (const auto r : buckets[(k*n[1] + j)*n[0] + i])
                        {
                            const Scalar d = (reference[r] - x).two_norm2();
                            if (d < best)
                            {
                                best = d;
                                nearest[p] = r;
                            }
                        }
                    }
            // all points closer than the searched rings have been seen
            if (best < std::numeric_limits<Scalar>::max() && ring*h >= std::sqrt(best))
                break;
            if (ring > n[0] + n[1] + n[2])
                break;
        }
    }
    return nearest;
}

/*!
 * \brief Pore pressures and fluid forces of a sphere packing in a box
 *
 * remesh() triangulates the packing with the walls and factorises the flow problem of the new network.
 * Between remeshes, update() takes the rates of change of the pore volumes from the motion of the spheres
 * and walls since the previous call, solves for the pore pressures and evaluates the pressure forces on
 * the spheres and walls (Chareyre et al. 2012; Catalano et al. 2014). For a compressible fluid each new
 * pore starts from the pressure of the nearest pore of the previous network.
 */
template<class Scalar>
class PoreScaleFiniteVolume
{
public:
    /*!
     * \param viscosity dynamic viscosity of the fluid
     * \param wallPressure imposed pressure of each wall (2*axis + side), none for an impermeable wall
     * \param bulkModulus bulk modulus of the fluid, none for an incompressible fluid
     * \param dt time step of all updates, needed for a compressible fluid
     */
    PoreScaleFiniteVolume(Scalar viscosity, const std::array<std::optional<Scalar>, 6>& wallPressure,
                          std::optional<Scalar> bulkModulus = std::nullopt, Scalar dt = 1.0)
    : viscosity_(viscosity), wallPressure_(wallPressure), bulkModulus_(bulkModulus), dt_(dt)
    {}

    void remesh(const std::vector<Point<Scalar>>& centers, const std::vector<Scalar>& radii, const Walls<Scalar>& walls)
    {
        std::vector<Point<Scalar>> oldPositions;
        std::vector<Scalar> oldPressure;
        if (network_ && bulkModulus_)
        {
            for (const auto& pore : network_->pores)
                oldPositions.push_back(pore.position);
            oldPressure = flow_->pressure();
        }

        radii_ = radii;
        triangulation_ = regularTriangulation(centers, radii, walls);
        network_ = std::make_unique<Network<Scalar>>(extractNetwork(triangulation_, centers, radii, walls));
        flow_ = std::make_unique<PoreScaleFlow<Scalar>>(*network_, viscosity_, wallPressure_, bulkModulus_, dt_);
        volumes_ = poreBulkVolumes(*network_, triangulation_, centers, radii, &walls);
        forces_.assign(centers.size() + 6, Point<Scalar>(0.0));
        ++numRemeshes_;

        if (!oldPositions.empty())
        {
            std::vector<Point<Scalar>> positions;
            for (const auto& pore : network_->pores)
                positions.push_back(pore.position);
            const auto nearest = nearestPoints(oldPositions, positions);
            std::vector<Scalar> pressure(positions.size());
            for (std::size_t i = 0; i < positions.size(); ++i)
                pressure[i] = oldPressure[nearest[i]];
            flow_->setPressure(pressure);
        }
    }

    /*!
     * \brief Pressures and forces for the motion of spheres and walls over the time step since the last call,
     *        with a fluid volume source per void volume (e.g. the expansion of freezing water)
     */
    void update(const std::vector<Point<Scalar>>& centers, const Walls<Scalar>& walls, Scalar dt, Scalar sourceRate = 0.0)
    {
        if (bulkModulus_ && dt != dt_)
            DUNE_THROW(Dune::InvalidStateException, "The time step of a compressible fluid is fixed at construction");
        const auto volumes = poreBulkVolumes(*network_, triangulation_, centers, radii_, &walls);
        std::vector<Scalar> rates(volumes.size()), sources(volumes.size());
        for (std::size_t i = 0; i < volumes.size(); ++i)
        {
            rates[i] = (volumes[i] - volumes_[i])/dt;
            sources[i] = sourceRate*network_->pores[i].volume;
        }
        volumes_ = volumes;
        flow_->solve(rates, sources);
        forces_ = fluidForces(*network_, flow_->pressure(), centers.size() + 6);
    }

    //! forces on the spheres, followed by the six walls
    const std::vector<Point<Scalar>>& forces() const { return forces_; }
    const std::vector<Scalar>& pressure() const { return flow_->pressure(); }
    const Network<Scalar>& network() const { return *network_; }
    std::size_t numRemeshes() const { return numRemeshes_; }

private:
    Scalar viscosity_;
    std::array<std::optional<Scalar>, 6> wallPressure_;
    std::optional<Scalar> bulkModulus_;
    Scalar dt_;
    std::vector<Scalar> radii_;
    Triangulation triangulation_;
    std::unique_ptr<Network<Scalar>> network_;
    std::unique_ptr<PoreScaleFlow<Scalar>> flow_;
    std::vector<Scalar> volumes_;
    std::vector<Point<Scalar>> forces_;
    std::size_t numRemeshes_ = 0;
};

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
