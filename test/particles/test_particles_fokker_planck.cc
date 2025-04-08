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

#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>

#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/initialize.hh>

#include <dumux/io/particles/writer.hh>
#include <dumux/particles/fokkerplanck.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux::Properties {

namespace TTag {

template<int dim>
struct ParticleAdvectionDiffusionTest
{
    using InheritsFrom = std::tuple<Tracer, CCTpfaModel>;
    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using SolutionDependentAdvection = std::false_type;
    using SolutionDependentMolecularDiffusion = std::false_type;
    using SolutionDependentHeatConduction = std::false_type;
    using Grid = Dune::YaspGrid<dim>;
};
} // end namespace TTag

} // end namespace Dumux::Properties

template<int dim>
void run()
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::ParticleAdvectionDiffusionTest<dim>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GlobalPosition = Dune::FieldVector<Scalar, Grid::dimensionworld>;

    GlobalPosition upperRight(1.0);
    std::array<int, Grid::dimension> cells;
    cells.fill(getParam<int>("Grid.Cells"));
    Dumux::GridManager<Grid> gridManager;
    gridManager.init(upperRight, cells);

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    std::vector<Scalar> density(gridGeometry->gridView().size(0), 0.0);
    std::vector<Scalar> densityAnalytical(gridGeometry->gridView().size(0), 0.0);

    const auto vel = Dune::FieldVector<Scalar, Grid::dimension>(getParam<Scalar>("FokkerPlanck.Velocity"));
    const auto center = Dune::FieldVector<Scalar, Grid::dimensionworld>(getParam<Scalar>("FokkerPlanck.DiracPosition"));
    const auto D = getParam<Scalar>("FokkerPlanck.DiffusionCoefficient");
    const auto updateAnalyticalSolution = [&](const Scalar t)
    {
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            const auto x = element.geometry().center() - center - vel*t;
            densityAnalytical[eIdx] =
                std::exp(-x*x/(4.0*D*t)) / (std::sqrt(std::pow(4.0*M_PI*D*t, Grid::dimension)));
        }
    };

    Dumux::IO::GridWriter gridWriter(
        IO::Format::pvd_with(IO::Format::vtp.with({.encoder = GridFormat::Encoding::raw})),
        gridGeometry->gridView(), "test_particles_fokker_planck_grid_" + std::to_string(dim) + "d");
    gridWriter.setCellField("c", [&](const auto& e){
        return density[gridGeometry->elementMapper().index(e)];
    }, GridFormat::Precision<Scalar>{});
    gridWriter.setCellField("c_exact", [&](const auto& e){
        return densityAnalytical[gridGeometry->elementMapper().index(e)];
    }, GridFormat::Precision<Scalar>{});
    gridWriter.write(0.0);

    // Run a comparison of the particle-based Fokker-Planck solver
    // with the analytical solution of the advection-diffusion equation

    // make a particle cloud
    using Particle = Dumux::Particles::Particle<Grid::dimensionworld, Scalar>;
    using Cloud = Dumux::Particles::SimpleParticleCloud<Particle>;
    auto cloud = std::make_shared<Cloud>();
    const auto numParticles = getParam<int>("Particles.NumParticles");
    cloud->resize(numParticles);

    Particles::FokkerPlanck<Scalar, GridGeometry> fokkerPlanck(gridGeometry, cloud);

    // initialize particles in one point (Dirac delta)
    const auto eIdx = intersectingEntityCartesianGrid(center, GlobalPosition(0.0), upperRight, cells);
    fokkerPlanck.init([&](auto& p, auto& particleData){
        p.activate();
        p.setPosition(center);
        particleData.eIdx = eIdx;
        particleData.time = 0.0;
    });

    // set constant velocity field
    std::vector<Dune::FieldVector<Scalar, Grid::dimension>> velocities(1, vel);
    fokkerPlanck.setVelocityData(velocities);

    Dumux::IO::Particles::Writer writer(fokkerPlanck.particleCloud(),
        "test_particles_fokker_planck_particles_" + std::to_string(dim) + "d"
    );

    const auto writeParticles = getParam<bool>("Output.WriteParticles", false);
    if (writeParticles)
        writer.write(0.0);

    const auto tEnd = getParam<double>("Particles.EndTime");
    const auto dt = getParam<double>("Particles.TimeStepSize");
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0, dt, tEnd);

    timeLoop->start();
    while (!timeLoop->finished())
    {
        fokkerPlanck.run(timeLoop->timeStepSize());

        timeLoop->advanceTimeStep();

        fokkerPlanck.computeDensity(density);
        updateAnalyticalSolution(timeLoop->time());

        Dune::Timer timer;
        if (writeParticles)
            writer.write(timeLoop->time());
        gridWriter.write(timeLoop->time());
        std::cout << "Writing took: " << timer.elapsed() << " seconds" << std::endl;

        timeLoop->reportTimeStep();
    }
    timeLoop->finalize();

    // print parameter report
    Parameters::print();
}

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init(argc, argv);

    const auto dim = Dumux::getParam<int>("Grid.Dimension");
    switch (dim)
    {
        case 1:
            run<1>();
            break;
        case 2:
            run<2>();
            break;
        case 3:
            run<3>();
            break;
        default:
            DUNE_THROW(Dune::IOError, "" << dim << "-dimensional case not implemented");
    }

    return 0;
} // end main
