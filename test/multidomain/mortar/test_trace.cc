#include <cmath>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/mortar/trace.hh>

#include "traceproblem.hh"
#include "tracespatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct Trace { using InheritsFrom = std::tuple<OneP>; };
struct TraceTpfa { using InheritsFrom = std::tuple<Trace, CCTpfaModel>; };
}  // namespace TTAG

template<class TypeTag>
struct Problem<TypeTag, TTag::Trace>
{ using type = TraceTestProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Trace>
{
    using type = TraceTestSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Scalar>
    >;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Trace>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> > ;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::Trace>
{ using type = Dune::YaspGrid<2>; };

}  // namespace Dumux::Properties

int main(int argc, char** argv) {
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    gridManager.loadBalance();
    const auto& gridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridView);
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x.resize(gridGeometry->numDofs());
    for (const auto& e : elements(gridGeometry->gridView()))
        x[gridGeometry->elementMapper().index(e)] = pressureField(e.geometry().center());
    gridVariables->init(x);

    FVTrace trace{gridVariables, [] (const auto&) { return true; }};
    using FV = GetPropType<TypeTag, Properties::FluxVariables>;
    const auto traceFluxes = trace.assembleFaceAverage(x, defaultCCAdvectiveFluxFunction<FV>(*problem));
    const auto tracePressures = trace.assembleFaceAverage(x, tracePressureFunctionTpfaOneP<FV>(*problem));

    // verify the result
    for (const auto& e : elements(trace.gridView())) {
        const auto center = e.geometry().center();
        const auto reconstructedPressure = tracePressures[trace.elementMapper().index(e)];
        const auto expectedPressure = pressureField(e.geometry().center());
        const auto pressureDiff = std::abs(expectedPressure - reconstructedPressure);
        if (pressureDiff > 1e-6)
            DUNE_THROW(Dune::InvalidStateException, "Pressure diff of " << pressureDiff << " (" << reconstructedPressure << " vs " << expectedPressure << ") at " << center << " detected.");

        const auto reconstructedFlux = traceFluxes[trace.elementMapper().index(e)];
        const auto exactGradient = pressureFieldGradient(e.geometry().center());
        const bool onHorizontalBoundary = center[1] < 1e-6 || center[1] > 1.0 - 1e-6;
        const auto expectedFlux = onHorizontalBoundary ? exactGradient[1] : exactGradient[0];
        const auto fluxDiff = std::abs(std::abs(expectedFlux) - std::abs(reconstructedFlux));
        if (fluxDiff > 1e-6)
            DUNE_THROW(Dune::InvalidStateException, "Flux diff of " << fluxDiff << " (" << reconstructedFlux << " vs " << expectedFlux << ") at " << center << " detected.");
    }

    return 0;
}
