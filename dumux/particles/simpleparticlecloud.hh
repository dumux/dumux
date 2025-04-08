// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Particles
 * \brief A simple implementation of a cloud of particles
 */
#ifndef DUMUX_PARTICLES_SIMPLE_PARTICLE_CLOUD_HH
#define DUMUX_PARTICLES_SIMPLE_PARTICLE_CLOUD_HH

#include <vector>
#include <ranges>

#include <dumux/parallel/parallel_for.hh>

#include <dumux/particles/particle.hh>

namespace Dumux::Particles {

/*!
 * \ingroup Particles
 * \brief A simple implementation of a cloud of particles
 *
 * The simple particle cloud assumes that particles are
 * independent of each other and can be updated in parallel.
 *
 * It allows to add and remove particles from the cloud
 * and modify particles in place.
 */
template<class P>
class SimpleParticleCloud
{
public:
    using Particle = P;

    SimpleParticleCloud() = default;

    /*!
     * \brief Resize the cloud to contain numParticles particles
     */
    void resize(std::size_t numParticles)
    {
        if (numParticles > particles_.size())
        {
            for (std::size_t pId = particles_.size(); pId < numParticles; ++pId)
                particles_.emplace_back(pId, false);
            if (particles_.size() != numParticles)
                DUNE_THROW(Dune::Exception, "");
        }
    }

    /*!
     * \brief the number of particles in the cloud
     */
    std::size_t size(bool onlyActive = true) const
    {
        if (onlyActive)
            return std::count_if(particles_.begin(), particles_.end(), [](const auto& p){ return p.isActive(); });
        else
            return particles_.size();
    }

    /*!
     * \brief Bring one inactive particle to life
     * Activates one particle from the pool of inactive particles
     * \note the particle is automatically activated
     * \note Take care to place it somewhere afterwards! Otherwise it's hanging somewhere in space.
     * \note If no particles are left new particles are added
     * \note This method is not thread safe
     */
    Particle& spawnParticle()
    {
        auto pIt = std::find_if(particles_.begin(), particles_.end(), [](const auto& p){ return !p.isActive(); });
        if (pIt == particles_.end())
        {
            std::cout << "-- running out of particles -- adding 10% more particles." << std::endl;
            const auto oldSize = particles_.size();
            this->resize(static_cast<std::size_t>(1.1*oldSize));
            pIt = particles_.begin() + oldSize;
        }

        pIt->activate();
        return *pIt;
    }

    /*!
     * \brief Deactivate a given particle
     */
    void deactivateParticle(const Particle& p)
    {
        particles_[p.id()].deactivate();
    }

    /*!
     * \brief Const range generator for iterating over all active particles
     * \note usage: for(const auto& particle : particles(cloud))
     * This is a free function found by means of ADL
     */
    friend inline std::ranges::forward_range auto particles(const SimpleParticleCloud& cloud)
    { return cloud.particles_ | std::views::filter([](const auto& p){ return p.isActive(); }); }

    /*!
     * \brief For each active particle apply a function
     * \note The function may modify the particle
     */
    template<class Func>
    void forEachActiveParticle(const Func& func)
    {
        parallelFor(particles_.size(), [&](std::size_t i){
            if (particles_[i].isActive())
                func(particles_[i]);
        });
    }

    /*!
     * \brief For each particle apply a function
     * \note The function may modify the particle
     */
    template<class Func>
    void forEachParticle(const Func& func)
    {
        parallelFor(particles_.size(), [&](std::size_t i){
            func(particles_[i]);
        });
    }
private:
    std::vector<Particle> particles_;
};

} // end namespace Dumux::Particles

#ifndef DOXYGEN // hide from doxygen
#ifdef DUMUX_HAVE_GRIDFORMAT
#include <gridformat/gridformat.hpp>

// make the SimpleParticleCloud compatible with GridFormat
namespace GridFormat::Traits {

template<typename P>
struct Cells<Dumux::Particles::SimpleParticleCloud<P>> {
    static std::ranges::forward_range auto get(const Dumux::Particles::SimpleParticleCloud<P>& c) {
        return particles(c);
    }
};

template<typename P>
struct Points<Dumux::Particles::SimpleParticleCloud<P>>{
    static std::ranges::forward_range auto get(const Dumux::Particles::SimpleParticleCloud<P>& c) {
        return particles(c);
    }
};

template<typename P>
struct CellPoints<Dumux::Particles::SimpleParticleCloud<P>, P> {
    static std::ranges::range auto get(const Dumux::Particles::SimpleParticleCloud<P>&, const P& p) {
        return std::views::single(p);
    }
};

template<typename P>
struct CellType<Dumux::Particles::SimpleParticleCloud<P>, P> {
    static auto get(const Dumux::Particles::SimpleParticleCloud<P>&, const P& p) {
        return GridFormat::CellType::vertex;
    }
};

template<typename P>
struct PointCoordinates<Dumux::Particles::SimpleParticleCloud<P>, P> {
    static GridFormat::Concepts::StaticallySizedRange auto get(const Dumux::Particles::SimpleParticleCloud<P>& c, const P& p) {
        return p.position();
    }
};

template<typename P>
struct PointId<Dumux::Particles::SimpleParticleCloud<P>, P> {
    static auto get(const Dumux::Particles::SimpleParticleCloud<P>& c, const P& p) {
        return p.id();
    }
};

template<typename P>
struct NumberOfPoints<Dumux::Particles::SimpleParticleCloud<P>> {
    static std::size_t get(const Dumux::Particles::SimpleParticleCloud<P>& c) {
        return c.size();
    }
};

template<typename P>
struct NumberOfCells<Dumux::Particles::SimpleParticleCloud<P>> {
    static std::size_t get(const Dumux::Particles::SimpleParticleCloud<P>& c) {
        return c.size();
    }
};

template<typename P>
struct NumberOfCellPoints<Dumux::Particles::SimpleParticleCloud<P>, P> {
    static std::size_t get(const Dumux::Particles::SimpleParticleCloud<P>& c, const P& p) {
        return 1;
    }
};

static_assert(GridFormat::Concepts::UnstructuredGrid<Dumux::Particles::SimpleParticleCloud<Dumux::Particles::Particle<1>>>);
static_assert(GridFormat::Concepts::UnstructuredGrid<Dumux::Particles::SimpleParticleCloud<Dumux::Particles::Particle<2>>>);
static_assert(GridFormat::Concepts::UnstructuredGrid<Dumux::Particles::SimpleParticleCloud<Dumux::Particles::Particle<3>>>);

} // end namespace GridFormat::Traits
#endif // DUMUX_HAVE_GRIDFORMAT
#endif // DOXYGEN

#endif
