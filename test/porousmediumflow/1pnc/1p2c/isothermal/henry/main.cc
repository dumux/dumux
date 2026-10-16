// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Henry problem: steady isochlors compared with Tables D1 and D2 of Fahs et al. \cite Fahs2016
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

//! Rows (Z, X10, X50, X90) of Table D1 or D2 of Fahs et al. \cite Fahs2016
std::vector<std::array<double, 4>> readIsochlorTable(const std::string& fileName)
{
    std::ifstream file(fileName);
    if (!file)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

    std::vector<std::array<double, 4>> rows;
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || !std::isdigit(static_cast<unsigned char>(line[0])))
            continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream stream(line);
        std::array<double, 4> row;
        for (auto& value : row)
            stream >> value;
        rows.push_back(row);
    }
    return rows;
}

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::HenryFahs;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize"));
    timeLoop->setPeriodicCheckPoint(getParam<Scalar>("TimeLoop.OutputInterval"));

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    VtkOutputModule vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        if (timeLoop->isCheckPoint() || timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // first crossing from the freshwater side on the grid line at height z
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    constexpr int soluteIdx = FluidSystem::comp1Idx;
    const auto& lowerLeft = gridGeometry->bBoxMin();
    const Scalar depth = gridGeometry->bBoxMax()[1] - lowerLeft[1];
    const auto isochlorPosition = [&](Scalar z, Scalar level)
    {
        std::vector<std::pair<Scalar, Scalar>> line;
        for (const auto& vertex : vertices(leafGridView))
        {
            const auto pos = vertex.geometry().corner(0);
            if (std::abs(pos[1] - z) < 1e-8*depth)
                line.emplace_back(pos[0], x[gridGeometry->vertexMapper().index(vertex)][soluteIdx]);
        }
        if (line.empty())
            DUNE_THROW(Dune::InvalidStateException, "No grid line at z = " << z);

        std::sort(line.begin(), line.end());
        for (std::size_t i = 1; i < line.size(); ++i)
        {
            const auto& [x0, c0] = line[i-1];
            const auto& [x1, c1] = line[i];
            if (c1 >= level)
                return x0 + (level - c0)*(x1 - x0)/(c1 - c0);
        }
        DUNE_THROW(Dune::InvalidStateException, "No isochlor " << level << " at z = " << z);
    };

    // the reference coordinates are in units of the domain depth
    const auto reference = readIsochlorTable(getParam<std::string>("Problem.ReferenceTable"));
    const std::array<Scalar, 3> levels = {0.1, 0.5, 0.9};
    Scalar maxError = 0.0;
    std::cout << std::fixed << std::setprecision(4) << "    Z      X10     X50     X90   (error in m)\n";
    for (const auto& row : reference)
    {
        std::cout << std::setw(7) << row[0];
        for (int i = 0; i < 3; ++i)
        {
            const Scalar position = (isochlorPosition(lowerLeft[1] + row[0]*depth, levels[i]) - lowerLeft[0])/depth;
            const Scalar error = (position - row[i+1])*depth;
            maxError = std::max(maxError, std::abs(error));
            std::cout << std::setw(8) << error;
        }
        std::cout << "\n";
    }
    std::cout << "Maximum isochlor error: " << maxError << " m" << std::endl;

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    if (!(maxError < getParam<Scalar>("Problem.IsochlorTolerance")))
        DUNE_THROW(Dune::Exception, "Maximum isochlor error " << maxError << " m exceeds the tolerance");

    return 0;
}
