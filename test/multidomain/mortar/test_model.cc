#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/multidomain/mortar/decomposition.hh>
#include <dumux/multidomain/mortar/model.hh>


template<typename MortarSolution, typename GridGeometry>
struct Solver : public Dumux::Mortar::Solver<MortarSolution, GridGeometry>
{
    using ParentType = Dumux::Mortar::Solver<MortarSolution, GridGeometry>;
    using ParentType::ParentType;

    void solve() override {}
};


int main(int argc, char** argv) {
    Dumux::initialize(argc, argv);

    using SubDomainGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
    using MortarGrid = Dune::FoamGrid<1, 2>;
    SubDomainGrid top{{0.0, 0.0}, {1.0, 1.0}, {10, 10}};
    SubDomainGrid bottom{{0.0, -1.0}, {1.0, 0.0}, {10, 10}};
    auto mortar = Dune::StructuredGridFactory<MortarGrid>::createCubeGrid({0.0, 0.0}, {1.0, 0.0}, {5});

    using SubDomainGridGeometry = Dumux::BoxFVGridGeometry<double, typename SubDomainGrid::LeafGridView>;
    using MortarGridGeometry = Dumux::BoxFVGridGeometry<double, typename MortarGrid::LeafGridView>;
    auto topGG = std::make_shared<SubDomainGridGeometry>(top.leafGridView());
    auto bottomGG = std::make_shared<SubDomainGridGeometry>(bottom.leafGridView());
    auto mortarGG = std::make_shared<MortarGridGeometry>(mortar->leafGridView());

    using MortarSolution = std::vector<double>;
    using SubDomainSolver = Solver<MortarSolution, SubDomainGridGeometry>;
    auto bottomSolver = std::make_shared<SubDomainSolver>(bottomGG);
    auto topSolver = std::make_shared<SubDomainSolver>(topGG);

    auto model = Dumux::Mortar::ModelFactory<MortarSolution, MortarGridGeometry, SubDomainGridGeometry>{}
        .withMortar(mortarGG)
        .withSubDomain(topSolver)
        .withSubDomain(bottomSolver)
        .make();

    int exitCode = 0;
    if (model.mortarSolution().size() != mortarGG->numDofs())
    {
        std::cout << "Unexpected mortar solution size" << std::endl;
        exitCode = 1;
    }

    return exitCode;
}
