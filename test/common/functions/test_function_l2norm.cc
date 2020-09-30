#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>
#include <dumux/multidomain/glue.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameters
    Dumux::Parameters::init([](auto& p){
        p["Grid.UpperRight"] = "1 1";
        p["Coarse.Grid.Cells"] = "10 10";
        p["Fine.Grid.Cells"] = "33 33";
    });

    // initialize a grid
    using Grid = Dune::YaspGrid<2>;
    Dumux::GridManager<Grid> gridManagerCoarse, gridManagerFine;
    gridManagerCoarse.init("Coarse");
    gridManagerFine.init("Fine");
    const auto gvFine = gridManagerFine.grid().leafGridView();
    const auto gvCoarse = gridManagerCoarse.grid().leafGridView();

    // get some types
    using R = Dune::FieldVector<double, 2>;
    using SolutionVector = Dune::BlockVector<R>;
    using GridView = Grid::LeafGridView;
    using namespace Dune::Functions;
    using namespace Dune::Functions::BasisFactory;
    using P2 = LagrangeBasis<GridView, 2>;

    // construct glue, any grid geometry should be fine here
    Dumux::BaseGridGeometry<GridView, Dumux::DefaultMapperTraits<GridView>> ggCoarse(gvCoarse), ggFine(gvFine);
    auto glue = makeGlue(ggFine, ggCoarse);

    // a quadratic function on the grid
    auto f = [&] (const auto& pos) -> R { return {pos.two_norm2(), pos[0]*pos[0]*pos[1]}; };

    // test Dumux::integrateGridFunction with unit function
    auto basisUnit = makeBasis(gvFine, lagrange<1>());

    using RUnit = Dune::FieldVector<double, 1>;
    Dune::BlockVector< RUnit > solUnit;
    interpolate(basisUnit, solUnit, [&] (const auto& pos) { return 1.0; });

    auto gfUnit = makeDiscreteGlobalBasisFunction<RUnit>(basisUnit, solUnit);
    auto integral = Dumux::integrateGridFunction(gvFine, gfUnit, 1);

    using std::abs;
    if ( abs(integral - 1.0) > 1e-13)
        DUNE_THROW(Dune::MathError, "Integral of grid function is not 1 (error = " << abs(integral - 1.0) << ")");

    // test Dumux::integrateL2Error
    SolutionVector solCoarse, solFine;
    interpolate(makeBasis(gvCoarse, power<2>(lagrange<2>())), solCoarse, f);
    interpolate(makeBasis(gvFine, power<2>(lagrange<2>())), solFine, f);

    P2 basisCoarse(gvCoarse), basisFine(gvFine);
    auto gfCoarse = makeDiscreteGlobalBasisFunction<R>(basisCoarse, solCoarse);
    auto gfFine = makeDiscreteGlobalBasisFunction<R>(basisFine, solFine);

    // project fine discrete function onto coarse grid and convert back to grid function
    const auto projector = Dumux::makeProjector(basisFine, basisCoarse, glue);
    auto params = projector.defaultParams();
    params.residualReduction = 1e-16;

    auto solCoarse2 = projector.project(solFine, params);
    auto gfCoarse2 = makeDiscreteGlobalBasisFunction<R>(basisCoarse, solCoarse2);

    // check that these functions are identical and that they exactly represent the analytical function
    auto l2norm = Dumux::integrateL2Error(gvCoarse, gfCoarse, gfCoarse2, 3);
    if (l2norm > 1e-14)
        DUNE_THROW(Dune::MathError, "L2 norm of projection and interpolation is non-zero (" << l2norm << ")");

    l2norm = Dumux::integrateL2Error(gvCoarse, makeAnalyticGridViewFunction(f, gvCoarse), gfCoarse2, 3);
    if (l2norm > 1e-14)
        DUNE_THROW(Dune::MathError, "L2 norm of analytical solution and interpolation is non-zero (" << l2norm << ")");

    // check the FV-style interface of integrateL2Error
    Dumux::BoxFVGridGeometry<double, GridView> ggBoxCoarse(gvCoarse), ggBoxFine(gvFine);
    auto f2 = [&] (const auto& pos) -> R { return {pos[0] + pos[1], pos[0]*pos[1]}; };

    interpolate(makeBasis(gvFine, power<2>(lagrange<1>())), solFine, f2);
    interpolate(makeBasis(gvCoarse, power<2>(lagrange<1>())), solCoarse, f2);

    // project fine discrete function onto coarse grid
    auto p1BasisCoarse = Dumux::getFunctionSpaceBasis(ggBoxCoarse);
    auto p1BasisFine = Dumux::getFunctionSpaceBasis(ggBoxFine);
    const auto projector2 = Dumux::makeProjector(p1BasisFine, p1BasisCoarse, glue);
    solCoarse2 = projector2.project(solFine, params);

    l2norm = Dumux::integrateL2Error(ggBoxCoarse, solCoarse, solCoarse2, 3);
    if (l2norm > 1e-14)
        DUNE_THROW(Dune::MathError, "L2 norm of FV-style projected solutions is non-zero (" << l2norm << ")");

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
