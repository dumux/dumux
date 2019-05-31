#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>
#include <dumux/multidomain/glue.hh>
#include <dumux/discretization/projection/projector.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

int main (int argc, char *argv[]) try
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
    using R = Dune::FieldVector<double, 1>;
    using SolutionVector = Dune::BlockVector<R>;
    using GridView = Grid::LeafGridView;
    using namespace Dune::Functions;
    using P2 = LagrangeBasis<GridView, 2>;

    // construct glue, any grid geometry should be fine here
    Dumux::BaseGridGeometry<GridView, Dumux::DefaultMapperTraits<GridView>> ggCoarse(gvCoarse), ggFine(gvFine);
    auto glue = makeGlue(ggFine, ggCoarse);

    // a quadratic function on the grid
    auto f = [](const auto& pos){ return pos.two_norm2(); };

    // make bases
    P2 basisCoarse(gvCoarse), basisFine(gvFine);
    SolutionVector solCoarse, solFine;
    interpolate(basisCoarse, solCoarse, f);
    interpolate(basisFine, solFine, f);
    auto gfCoarse = makeDiscreteGlobalBasisFunction<R>(basisCoarse, solCoarse);
    auto gfFine = makeDiscreteGlobalBasisFunction<R>(basisFine, solFine);

    // project fine discrete function onto coarse grid and convert back to grid function
    const auto projector = Dumux::makeProjector(basisFine, basisCoarse, glue);
    SolutionVector solCoarse2;
    projector.project(solFine, solCoarse2);
    auto gfCoarse2 = makeDiscreteGlobalBasisFunction<R>(basisCoarse, solCoarse);

    // check that these functions are identical and that they exactly represent the analytical function
    auto l2norm = Dumux::integrateL2Error(gvCoarse, gfCoarse, gfCoarse2, 3);
    if (l2norm > 1e-14)
        DUNE_THROW(Dune::MathError, "L2 norm of projection and interpolation is non-zero (" << l2norm << ")");

    l2norm = Dumux::integrateL2Error(gvCoarse, makeAnalyticGridViewFunction(f, gvCoarse), gfCoarse2, 3);
    if (l2norm > 1e-14)
        DUNE_THROW(Dune::MathError, "L2 norm of analytical solution and interpolation is non-zero (" << l2norm << ")");

    // check the FV-style interface of integrateL2Error
    Dumux::BoxFVGridGeometry<double, GridView> ggBoxCoarse(gvCoarse), ggBoxFine(gvFine);
    auto p1BasisCoarse = Dumux::getFunctionSpaceBasis(ggBoxCoarse);
    auto p1BasisFine = Dumux::getFunctionSpaceBasis(ggBoxFine);
    auto f2 = [](const auto& pos){ return pos[0]; };
    interpolate(p1BasisFine, solFine, f2);
    interpolate(p1BasisCoarse, solCoarse, f2);

    // project fine discrete function onto coarse grid
    const auto projector2 = Dumux::makeProjector(p1BasisFine, p1BasisCoarse, glue);
    projector2.project(solFine, solCoarse2);

    l2norm = Dumux::integrateL2Error(ggBoxCoarse, solCoarse, solCoarse2, 3);
    if (l2norm > 1e-14)
        DUNE_THROW(Dune::MathError, "L2 norm of FV-style projected solutions is non-zero (" << l2norm << ")");

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
