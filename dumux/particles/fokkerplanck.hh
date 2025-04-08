// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Particles
 * \copydoc Dumux::Particles::FokkerPlanck
 */
#include <memory>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>


#include <dumux/parallel/parallel_for.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/geometryintersection.hh>

#include <dumux/particles/particle.hh>
#include <dumux/particles/simpleparticlecloud.hh>

namespace Dumux::Particles {

/*!
 * \ingroup Particles
 * \ingroup FokkerPlanck
 *
 * \brief A particle solver for the Fokker-Planck equation
 *
 * A cloud of particles, each particle associated with a given mass \f$m_i\f$,
 * is evolved such that the probability density function \f$\rho(\mathbf{x}, t)\f$
 * of the particle position \f$ \mathbf{X}_t \f$
 * approximates the solution of the Fokker-Planck equation:
 *
 * \f[
 * \frac{\partial \rho}{\partial t} + \nabla \cdot (\boldsymbol{\mu} \rho + D \nabla \rho) = 0
 * \f]
 *
 * Particles are evolved with the \ref Dumux::Particles::FokkerPlanck::run method,
 * which updates each particles position \f$\mathbf{X}_t\f$ with the Euler scheme:
 *
 * \f[
 * \Delta\mathbf{X}_t = \boldsymbol{\mu}(\mathbf{X}_t, t) \Delta t + \sqrt{2 D} \Delta \mathbf{W}_t
 * \f]
 *
 * The particle density on the grid can be computed with
 * the \ref Dumux::Particles::FokkerPlanck::computeDensity method as
 *
 * \f[
 * \rho_E = \frac{1}{|E|} \sum_{i=1}^{N_E} m_i,
 * \f]
 *
 * where \f$|E|\f$ is the volume of the grid element \f$E\f$ and
 * \f$m_i\f$ is the mass of the particle \f$i\f$ in element \f$E\f$.
 *
 * The Fokker-Planck equation for particle positions is equivalent to the
 * advection-diffusion equation for the particle density.
 *
 *
 * \todo the current policy when particles cross the grid boundary
 *       is to deactivate them. We may want to add a customizable
 *       boundary policy in the future.
 *
 * \todo so far, the only supports isotropic
 *       and spatially-independent diffusion tensors
 *
 * \todo we currently assume the velocity field in given as a
 *       piecewise-constant function on the grid elements. This should
 *       be generalized to a piecewise polynomials function in the future.
 */
template<class Scalar, class GridGeometry>
class FokkerPlanck
{
    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using Particle = Dumux::Particles::Particle<dimWorld, Scalar>;
    using Cloud = SimpleParticleCloud<Particle>;
    using GlobalPosition = Particle::GlobalPosition;
    using Velocity = Dune::FieldVector<Scalar, dim>;

    static_assert(GridGeometry::discMethod == DiscretizationMethods::cctpfa, "Implementation currently assumes a CC-TPFA grid geometry");
public:
    FokkerPlanck(std::shared_ptr<const GridGeometry> gg, std::shared_ptr<Cloud> cloud)
    : gridGeometry_(std::move(gg))
    , cloud_(std::move(cloud))
    , rndGen_(std::random_device{}())
    {
        if (cloud_->size(false) == 0)
            cloud_->resize(getParam<int>("FokkerPlanck.NumParticles"));

        const auto numParticles = cloud_->size(false);
        diffusionCoefficient_ = getParam<Scalar>("FokkerPlanck.DiffusionCoefficient");
        data_.resize(numParticles);
        particleMass_ = getParam<Scalar>("FokkerPlanck.ParticleMass", 1.0/numParticles);
        volumes_.resize(gridGeometry_->gridView().size(0), 0.0);
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            const auto eIdx = gridGeometry_->elementMapper().index(element);
            const auto geometry = element.geometry();
            volumes_[eIdx] = geometry.volume();
        }
    }

    /*!
     * \brief data associated with each particle id
     * \note This is used to attach a particle to a grid element
     *       when initializing the particle cloud (see \ref Dumux::Particles::FokkerPlanck::init)
     */
    struct ParticleData
    {
        std::size_t eIdx; // in which element of the grid this particle is
        Scalar time; // at which time the particle lives
    };

    /*!
     * \brief For each particle apply a function to initialize it's state
     * \note The function may modify the particle and its associated data
     */
    template<class Func>
    void init(const Func& func)
    {
        cloud_->forEachParticle([&](auto& p){
            func(p, data_[p.id()]);
        });
    }

    /*!
     * \brief Evolve particle state with Euler scheme
     * \param dt the time step size
     *
     * Updates the particle position \f$\mathbf{X}_t\f$ by
     *
     * \f[
     * \Delta\mathbf{X}_t = \boldsymbol{\mu}(\mathbf{X}_t, t) \Delta t + \sqrt{2 D} \Delta \mathbf{W}_t
     * \f]
     *
     * where \f$\boldsymbol{\mu}(\mathbf{X}_t, t)\f$ is the velocity field,
     * \f$D\f$ is the diffusion coefficient, and \f$\Delta \mathbf{W}_t\f$
     * is a Wiener increment. Note that \f$\Delta \mathbf{W}_t\f$ is normal-distributed
     * with mean zero and variance \f$\Delta t\f$: \f$\Delta \mathbf{W}_t \sim \mathcal{N}(0, \Delta t)\f$.
     */
    void run(const Scalar dt)
    {
        if (velocities_.empty())
            DUNE_THROW(Dune::InvalidStateException, "A velocity field must be set before calling run");

        // move the particles with the given velocity field through the grid
        cloud_->forEachActiveParticle([&](auto& p){
            moveParticle_(p, dt);
        });
    }

    /*!
     * \brief computes the particle density in the grid
     * \param density the density data to fill for each grid element
     *
     * Computes a discrete piece-wise constant approximation of the density
     * with cell values given as
     *
     * \f[
     * \rho_E = \frac{1}{|E|} \sum_{i=1}^{N_E} m_i,
     * \f]
     *
     * where \f$|E|\f$ is the volume of the grid element \f$E\f$ and
     * \f$m_i\f$ is the mass of the particle \f$i\f$ in element \f$E\f$.
     */
    void computeDensity(std::vector<Scalar>& density) const
    {
        density.assign(gridGeometry_->gridView().size(0), 0.0);
        for (const auto& p : particles(*cloud_))
        {
            const auto eIdx = data_[p.id()].eIdx;
            density[eIdx] += particleMass_/volumes_[eIdx];
        }
    }

    /*!
     * \brief set the velocity data for each grid element
     * \param velocities the velocity data
     * \note to set a constant velocity field, use a vector of size 1
     */
    void setVelocityData(const std::vector<Velocity>& velocities)
    {
        if (velocities.size() != gridGeometry_->gridView().size(0) && velocities.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Velocity data size does not match number of elements or 1");
        velocities_ = velocities;
    }

    //! access the underlying particle cloud
    std::shared_ptr<const Cloud> particleCloud() const
    { return cloud_; }

private:
    void moveParticle_(Particle& p, const Scalar dt)
    {
        auto& particleData = data_[p.id()];

        // (1) Drift: time loop until time integration is finished
        particleData.time = dt; // recharge time
        int count = 0;
        while (particleData.time > 0.0 && velocity_(particleData.eIdx).two_norm2() > 0.0)
        {
            const auto tStartIteration = particleData.time;

            const auto velocityInGrid = velocity_(particleData.eIdx);
            moveParticleInGridElement_(p, velocityInGrid);

            count = particleData.time < tStartIteration ? 0 : count + 1;
            if (count > 4)
            {
                std::cout << "Particle [" << p.id() << "] at " << p.position() << " in element " << particleData.eIdx << ", t=" << particleData.time << ", velocity: " << velocityInGrid << "\n";
                std::cout << "-- during advection step: particle got stuck and is deactivated\n";

                cloud_->deactivateParticle(p);
                particleData.time = -1.0;
            }
        }

        static_assert(dim == dimWorld, "Not implemented yet for dim != dimWorld");

        // (2) Brownian motion: time loop until time integration is finished
        particleData.time = 1.0; // recharge time
        const auto velocityInGrid = [&]{
            const auto scale = std::sqrt(2.0*diffusionCoefficient_);
            Dune::FieldVector<Scalar, dimWorld> newPos(p.position());
            for (int i = 0; i < dim; ++i)
                newPos[i] += scale * std::normal_distribution<Scalar>(0.0, std::sqrt(dt))(rndGen_);
            return newPos - p.position();
        }();

        const auto velocityNorm = velocityInGrid.two_norm2();
        count = 0;
        while (particleData.time > 0.0 && velocityNorm > 0.0)
        {
            const auto tStartIteration = particleData.time;

            moveParticleInGridElement_(p, velocityInGrid);

            count = particleData.time < tStartIteration ? 0 : count + 1;
            if (count > 4)
            {
                std::cout << "Particle [" << p.id() << "] at " << p.position() << " in element " << particleData.eIdx << ", t=" << particleData.time << ", velocity: " << velocityInGrid << "\n";
                std::cout << "-- during diffusion step: particle got stuck and is deactivated\n";

                cloud_->deactivateParticle(p);
                particleData.time = -1.0;
            }
        }
    }

    void moveParticleInGridElement_(Particle& p, const Velocity& velocity)
    {
        auto& particleData = data_[p.id()];
        assert(particleData.time > 0.0);

        const auto eIdx = particleData.eIdx;
        if (eIdx >= gridGeometry_->gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "Element index is invalid");

        const auto& element = gridGeometry_->element(eIdx);
        const auto geometry = element.geometry();

        const auto startPos = p.position();
        auto endPos = startPos;
        endPos.axpy(particleData.time, velocity);

        // we stay inside this element
        if (intersectsPointGeometry(endPos, geometry))
        {
            p.setPosition(endPos);
            particleData.time = -1.0; // this particle is done
            return;
        }

        // increase the change of finding the intersection
        // if we are right on the element facet by backtracking a bit
        auto beforeStart = startPos;
        beforeStart.axpy(-1e-4*particleData.time, velocity);

        const auto fvGeometry = localView(*gridGeometry_).bind(element);
        Dune::AffineGeometry<double, 1, dimWorld> segGeo(Dune::GeometryTypes::line, std::array<GlobalPosition, 2>{{beforeStart, endPos}});
        using IP = Dumux::IntersectionPolicy::PointPolicy<double, dimWorld>;
        using IntersectionAlgorithm = GeometryIntersection<Dune::AffineGeometry<double, 1, dimWorld>, decltype(fvGeometry.geometry(fvGeometry.scvf(0))), IP>;
        using Intersection = typename IntersectionAlgorithm::Intersection;
        Intersection intersection;

        static_assert(dim >= 1, "Particle transport only makes sense for 1D and higher");
        const auto intersectionFound = [&](const auto& scvf, Intersection& intersection)
        {
            if constexpr (dim == 1)
            {
                const bool found = intersectsPointGeometry(scvf.ipGlobal(), segGeo);
                if (found) intersection = scvf.ipGlobal();
                return found;
            }
            else
                return IntersectionAlgorithm::intersection(
                    segGeo, fvGeometry.geometry(scvf), intersection
                );
        };

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (intersectionFound(scvf, intersection))
            {
                // for now we deactivate the particle if it hits a boundary
                if (scvf.boundary())
                {
                    cloud_->deactivateParticle(p);
                    particleData.time = -1.0; // this particle is done
                    return;
                }

                // corner case (intersects more than one scvfs)
                if (scvf.unitOuterNormal()*velocity < 1e-4*velocity.two_norm())
                    continue;

                // put particle on face into the neighboring element
                const auto velocityNorm = velocity.two_norm2();
                const auto distance = (intersection-startPos).two_norm2();
                intersection.axpy(1e-3, endPos-startPos);
                p.setPosition(intersection);
                particleData.eIdx = scvf.outsideScvIdx();
                particleData.time -= std::sqrt(distance/velocityNorm);
                return;
            }
        }

        // particle got lost: try to find the particle in the grid
        static constexpr bool isCartesianGrid = Dune::Capabilities::isCartesian<typename GridGeometry::GridView::Grid>::v;
        const auto entities = intersectingEntities(p.position(), gridGeometry_->boundingBoxTree(), isCartesianGrid);
        if (entities.empty()) // outside the grid
        {
            cloud_->deactivateParticle(p);
            particleData.time = -1.0; // this particle is done
            return;
        }

        // found particle in an element, try to move on from there
        particleData.eIdx = entities[0];
    }

    const Velocity& velocity_(std::size_t eIdx) const
    {
        if (velocities_.size() == 1)
            return velocities_[0];
        return velocities_[eIdx];
    }

    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::vector<ParticleData> data_;
    std::shared_ptr<Cloud> cloud_;
    std::mt19937 rndGen_;

    std::vector<Velocity> velocities_;
    Scalar diffusionCoefficient_;

    std::vector<Scalar> volumes_; //! element volumes
    Scalar particleMass_;

};

} // end namespace Dumux::Particles
