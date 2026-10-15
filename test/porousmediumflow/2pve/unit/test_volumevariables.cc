// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
#include <config.h>

#include <array>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

#include "../cartesianinjection/properties.hh"

namespace {

template<class V>
concept HasSolidState = requires(const V& v) { v.solidState(); };

void checkClose(double actual, double expected, const char* quantity)
{
    if (std::abs(actual - expected) > 1e-12*std::abs(expected))
        DUNE_THROW(Dune::Exception, quantity << ": expected " << expected << ", got " << actual);
}

// Height scales integrated storage and permeability while averages and temperature stay fixed.
void checkColumn(double height)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TwoPVEImmiscibleTpfa;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FineGrid = GetPropType<TypeTag, Properties::FineGrid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FineGridGeometry = GetPropType<TypeTag, Properties::FineGridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    static_assert(!HasSolidState<VolumeVariables>);

    Grid coarseGrid(Dune::FieldVector<double, 1>({0.0}), Dune::FieldVector<double, 1>({4.0}), std::array<int, 1>({2}));
    FineGrid fineGrid(Dune::FieldVector<double, 2>({0.0, 0.0}), Dune::FieldVector<double, 2>({4.0, height}), std::array<int, 2>({2, 2}));
    auto coarseGeometry = std::make_shared<GridGeometry>(coarseGrid.leafGridView());
    auto fineGeometry = std::make_shared<FineGridGeometry>(fineGrid.leafGridView());
    auto fineProblem = std::make_shared<typename Problem::FineProblem>(fineGeometry);
    auto view = std::make_shared<typename Problem::FineLevelView>(fineGeometry, coarseGeometry, fineProblem->spatialParamsPtr());
    auto problem = std::make_shared<Problem>(coarseGeometry, view, fineProblem);
    GetPropType<TypeTag, Properties::SolutionVector> solution(coarseGeometry->numDofs());
    for (auto& priVars : solution)
    {
        priVars[0] = 1e7;
        priVars[1] = 0.1;
    }
    GetPropType<TypeTag, Properties::LocalResidual> residual(problem.get());
    for (const auto& element : elements(coarseGeometry->gridView()))
    {
        const auto fvGeometry = localView(*coarseGeometry).bindElement(element);
        const auto elemSol = elementSolution(element, solution, *coarseGeometry);
        for (const auto& scv : scvs(fvGeometry))
        {
            VolumeVariables volVars;
            volVars.update(elemSol, *problem, element, scv);
            checkClose(volVars.porosity(), 0.2, "average porosity");
            checkClose(volVars.porosityIntegral(), 0.2*height, "porosity integral");
            checkClose(volVars.extrusionFactor(), height, "extrusion height");
            checkClose(volVars.permeability(), 2e-12, "average permeability");
            checkClose(volVars.permeabilityIntegral(), 2e-12*height, "permeability integral");
            checkClose(volVars.temperature(), 326.0, "prescribed temperature");
            const auto storage = residual.computeStorage(*problem, scv, volVars);
            for (int phaseIdx = 0; phaseIdx < 2; ++phaseIdx)
            {
                checkClose(volVars.fluidState().temperature(phaseIdx), 326.0, "phase temperature");
                checkClose(storage[phaseIdx]*volVars.extrusionFactor(), 0.2*height*volVars.density(phaseIdx)*volVars.saturation(phaseIdx), "integrated phase storage");
                if (!std::isfinite(storage[phaseIdx]) || storage[phaseIdx] <= 0.0)
                    DUNE_THROW(Dune::Exception, "Phase storage must be finite and positive");
            }
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init([](auto& params) {
        params["Problem.Name"] = "test_2pve_volumevariables";
        params["Problem.EnableGravity"] = "true";
        params["BoundaryConditions.InjectionRate"] = "-1e-4";
        params["SpatialParams.Permeability"] = "2e-12";
        params["SpatialParams.Porosity"] = "0.2";
        params["SpatialParams.EnableLens"] = "false";
        params["SpatialParams.Swr"] = "0";
        params["SpatialParams.Snr"] = "0";
        params["SpatialParams.BrooksCoreyPcEntry"] = "1e5";
        params["SpatialParams.BrooksCoreyLambda"] = "2";
    });
    checkColumn(8.0); // The porosity integral exceeds one, while porosity remains dimensionless.
    checkColumn(16.0);
}
