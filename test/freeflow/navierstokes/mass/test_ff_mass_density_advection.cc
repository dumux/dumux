// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the upwinded density in the mass balance of the Navier-Stokes model.
 *
 * A density step of an ideal gas is advected with a prescribed uniform velocity, once in positive and once
 * in negative x-direction, the second setup being the mirror image of the first. Implicit Euler with a fully
 * upwinded density is monotone and does not depend on the flow direction: both solutions stay within the
 * initial density range and are mirror images of each other. The error with respect to the advected step
 * decreases with the grid size.
 */

#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/cvfelocalresidual_.hh>
#include <dumux/discretization/box.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

//! Density step moving with the flow, starting at the step position measured from the inflow boundary
template<class Scalar>
class DensityStep
{
public:
    explicit DensityStep(Scalar velocity)
    : velocity_(velocity)
    {
        length_ = getParam<Scalar>("Grid.Length");
        stepPosition_ = getParam<Scalar>("Problem.StepPosition");
        stepWidth_ = getParam<Scalar>("Problem.StepWidth");
        upstreamDensity_ = getParam<Scalar>("Problem.UpstreamDensity");
        downstreamDensity_ = getParam<Scalar>("Problem.DownstreamDensity");
    }

    Scalar density(Scalar x, Scalar time) const
    {
        using std::abs;
        using std::tanh;
        const Scalar distanceFromInflow = velocity_ > 0.0 ? x : length_ - x;
        const Scalar downstreamDistance = distanceFromInflow - stepPosition_ - abs(velocity_)*time;
        return upstreamDensity_ + 0.5*(downstreamDensity_ - upstreamDensity_)*(1.0 + tanh(downstreamDistance/stepWidth_));
    }

    Scalar velocity() const
    { return velocity_; }

    Scalar upstreamDensity() const
    { return upstreamDensity_; }

private:
    Scalar velocity_;
    Scalar length_;
    Scalar stepPosition_;
    Scalar stepWidth_;
    Scalar upstreamDensity_;
    Scalar downstreamDensity_;
};

template<class TypeTag>
class DensityStepProblem : public CVFENavierStokesMassProblem<TypeTag>
{
    using ParentType = CVFENavierStokesMassProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using typename ParentType::BoundaryTypes;
    using typename ParentType::BoundaryFluxes;

    DensityStepProblem(std::shared_ptr<const GridGeometry> gridGeometry, const DensityStep<Scalar>& step)
    : ParentType(gridGeometry)
    , step_(step)
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    //! The upstream density enters at the inflow boundary, the density of the solution leaves at the outflow boundary
    template<class ElementVariables, class FaceIpData>
    BoundaryFluxes boundaryFlux(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const FaceIpData& faceIpData) const
    {
        const Scalar normalVelocity = velocityAtPos(faceIpData.global())*faceIpData.unitOuterNormal();
        Scalar density = step_.upstreamDensity();
        if (normalVelocity > 0.0)
        {
            const auto& shapeValues = cache(elemVars, faceIpData).shapeValues();
            density = 0.0;
            for (const auto& localDof : localDofs(elemDisc))
                density += elemVars[localDof.index()].density()*shapeValues[localDof.index()][0];
        }

        BoundaryFluxes values(0.0);
        values[0] = density*normalVelocity;
        return values;
    }

    GlobalPosition velocityAtPos(const GlobalPosition&) const
    { return {step_.velocity(), 0.0}; }

private:
    DensityStep<Scalar> step_;
};

namespace Properties {

namespace TTag {
struct DensityStepBox { using InheritsFrom = std::tuple<NavierStokesMassOneP, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::DensityStepBox>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DensityStepBox>
{ using type = DensityStepProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DensityStepBox>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar>>;
};

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::DensityStepBox> { static constexpr bool value = true; };

} // end namespace Properties

using TypeTag = Properties::TTag::DensityStepBox;
using Scalar = GetPropType<TypeTag, Properties::Scalar>;

static_assert(std::is_base_of_v<Experimental::CVFELocalResidual<TypeTag>, GetPropType<TypeTag, Properties::LocalResidual>>,
              "The test has to assemble the flux integral of the variables per local dof");

//! Position and density of each degree of freedom at the end time
std::vector<std::array<Scalar, 3>> advectDensityStep(Scalar velocity, int numCellsX)
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Air = Components::Air<Scalar>;

    const auto length = getParam<Scalar>("Grid.Length");
    const auto height = getParam<Scalar>("Grid.Height");
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto cfl = getParam<Scalar>("TimeLoop.Cfl");
    const auto temperature = getParam<Scalar>("SpatialParams.Temperature");

    const Grid grid({length, height}, {numCellsX, 1});
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    const DensityStep<Scalar> step(velocity);
    auto problem = std::make_shared<Problem>(gridGeometry, step);

    SolutionVector x(gridGeometry->numDofs());
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        elemDisc.bind(element);
        for (const auto& scv : scvs(elemDisc))
            x[scv.dofIndex()] = Air::gasPressure(temperature, step.density(scv.dofPosition()[0], 0.0));
    }
    auto xOld = x;

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using std::abs;
    using std::ceil;
    const int numSteps = static_cast<int>(ceil(tEnd*abs(velocity)*numCellsX/(cfl*length)));
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, tEnd/numSteps, tEnd, false);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, std::make_shared<LinearSolver>());
    nonLinearSolver.setVerbosity(0);

    timeLoop->start(); do
    {
        nonLinearSolver.solve(x);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
    } while (!timeLoop->finished());

    std::vector<std::array<Scalar, 3>> densities(gridGeometry->numDofs());
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        elemDisc.bind(element);
        for (const auto& scv : scvs(elemDisc))
        {
            const auto& pos = scv.dofPosition();
            densities[scv.dofIndex()] = {pos[0], pos[1], Air::gasDensity(temperature, x[scv.dofIndex()][0])};
        }
    }
    return densities;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init(argc, argv, [](Dune::ParameterTree& params)
    {
        params["Grid.Length"] = "1.0";
        params["Grid.Height"] = "0.1";
        params["TimeLoop.TEnd"] = "0.3";
        params["TimeLoop.Cfl"] = "0.5";
        params["Problem.Velocity"] = "1.0";
        params["Problem.StepPosition"] = "0.3";
        params["Problem.StepWidth"] = "0.05";
        params["Problem.UpstreamDensity"] = "10.0";
        params["Problem.DownstreamDensity"] = "1.0";
        params["Problem.EnableGravity"] = "false";
        params["SpatialParams.Temperature"] = "293.15";
        params["Study.BaseCells"] = "20";
        params["Study.NumLevels"] = "3";
        params["Study.MaxMirrorDifference"] = "1e-10";
        params["Study.MaxOvershoot"] = "1e-10";
        params["Study.MinOrder"] = "0.6";
    });

    const auto length = getParam<Scalar>("Grid.Length");
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto velocity = getParam<Scalar>("Problem.Velocity");
    const auto upstreamDensity = getParam<Scalar>("Problem.UpstreamDensity");
    const auto downstreamDensity = getParam<Scalar>("Problem.DownstreamDensity");
    const auto baseCells = getParam<int>("Study.BaseCells");
    const auto numLevels = getParam<int>("Study.NumLevels");
    const auto maxMirrorDifference = getParam<Scalar>("Study.MaxMirrorDifference");
    const auto maxOvershoot = getParam<Scalar>("Study.MaxOvershoot");
    const auto minOrder = getParam<Scalar>("Study.MinOrder");

    using std::abs;
    const Scalar densityJump = abs(upstreamDensity - downstreamDensity);
    const Scalar minDensity = std::min(upstreamDensity, downstreamDensity);
    const Scalar maxDensity = std::max(upstreamDensity, downstreamDensity);

    const DensityStep<Scalar> exactStep(velocity);

    bool passed = true;
    Scalar previousError = 0.0;
    std::cout << std::setprecision(6) << std::scientific;
    for (int level = 0; level < numLevels; ++level)
    {
        const int numCellsX = baseCells << level;
        auto forward = advectDensityStep(velocity, numCellsX);
        auto backward = advectDensityStep(-velocity, numCellsX);

        // mirror the backward flow at the center of the domain
        for (auto& dof : backward)
            dof[0] = length - dof[0];

        const auto byPosition = [](const auto& a, const auto& b)
        { return a[0] < b[0] - 1e-12 || (abs(a[0] - b[0]) <= 1e-12 && a[1] < b[1]); };
        std::sort(forward.begin(), forward.end(), byPosition);
        std::sort(backward.begin(), backward.end(), byPosition);

        Scalar mirrorDifference = 0.0;
        Scalar overshoot = 0.0;
        Scalar error = 0.0;
        for (std::size_t i = 0; i < forward.size(); ++i)
        {
            const Scalar exact = exactStep.density(forward[i][0], tEnd);
            mirrorDifference = std::max(mirrorDifference, abs(forward[i][2] - backward[i][2])/densityJump);
            for (const Scalar density : {forward[i][2], backward[i][2]})
                overshoot = std::max({overshoot, (density - maxDensity)/densityJump, (minDensity - density)/densityJump});
            error += abs(forward[i][2] - exact)/densityJump;
        }
        error /= forward.size();

        std::cout << "[Result] cells " << numCellsX
                  << " mirror difference " << mirrorDifference
                  << " overshoot " << overshoot
                  << " mean error " << error;
        if (level > 0)
            std::cout << " order " << std::log2(previousError/error);
        std::cout << std::endl;

        if (!(mirrorDifference <= maxMirrorDifference))
        {
            std::cout << "[Failed] the solutions for both flow directions are not mirror images" << std::endl;
            passed = false;
        }
        if (!(overshoot <= maxOvershoot))
        {
            std::cout << "[Failed] the density leaves the initial range" << std::endl;
            passed = false;
        }
        if (level == numLevels - 1 && level > 0 && !(std::log2(previousError/error) >= minOrder))
        {
            std::cout << "[Failed] observed order is below " << minOrder << std::endl;
            passed = false;
        }
        previousError = error;
    }

    return passed ? 0 : 1;
}
