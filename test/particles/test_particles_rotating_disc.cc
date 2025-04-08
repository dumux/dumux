// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <memory>
#include <cmath>
#include <random>
#include <iostream>

#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/initialize.hh>

#include <dumux/particles/particle.hh>
#include <dumux/particles/simpleparticlecloud.hh>
#include <dumux/io/particles/writer.hh>

//! initialize particles on a disc
template<class Cloud>
void initialize(Cloud& cloud)
{
    const double radius = 4.0;
    const std::size_t numParticles = 100;
    cloud.resize(numParticles);

    std::mt19937 gen(0.0);
    std::uniform_real_distribution<double> distRadiusSq(0.0, radius*radius);
    std::uniform_real_distribution<double> distAngle(0.0, 2.0*M_PI);
    for (std::size_t pIdx = 0; pIdx < numParticles; ++pIdx)
    {
        auto& p = cloud.spawnParticle();
        const auto radius = std::sqrt(distRadiusSq(gen));
        const auto angle = distAngle(gen);
        p.setPosition({radius*std::cos(angle), radius*std::sin(angle)});
    }
}

//! evolve particles on rotating disc
template<class Cloud>
void evolve(Cloud& cloud, const double dt)
{
    const double angularVelocity = 2.0*M_PI/5.0;
    std::mt19937 gen(0.0);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    cloud.forEachActiveParticle([&](auto& p)
    {
        const auto& pos = p.position();
        const auto radius = std::hypot(pos[0], pos[1]) + 0.1*dist(gen);
        const auto angle = std::atan2(pos[1], pos[0]) + dt*angularVelocity;
        p.setPosition({radius*std::cos(angle), radius*std::sin(angle)});

        // deactivate the particle with a probability of 1%
        if (dist(gen) < 0.05)
            cloud.deactivateParticle(p);
    });
}

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    static constexpr int dimWorld = 2;
    using Particle = Dumux::Particles::Particle<dimWorld>;
    using Cloud = Dumux::Particles::SimpleParticleCloud<Particle>;
    auto cloud = std::make_shared<Cloud>();
    initialize(*cloud);

    // vtk output
    Dumux::IO::Particles::Writer writer(cloud, "test_particles_rotating_disc");
    std::vector<double> radius = Dumux::linspace(0.5, 0.7, cloud->size());
    writer.addParticleData(radius, "radius");
    writer.write(0.0);

    Dumux::TimeLoop<double> timeLoop(0.0, 0.1, 5.0);
    timeLoop.start();
    while (!timeLoop.finished())
    {
        evolve(*cloud, timeLoop.timeStepSize());
        timeLoop.advanceTimeStep();
        writer.write(timeLoop.time());
        timeLoop.reportTimeStep();

        std::cout << "-- number of active particles: " << cloud->size() << "\n"
                  << "-- number of inactive particles: " << cloud->size(false)-cloud->size(true) << std::endl;
    }

    timeLoop.finalize();

    return 0;
}
