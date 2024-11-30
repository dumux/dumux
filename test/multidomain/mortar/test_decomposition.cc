#include <numeric>

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/multidomain/mortar/decomposition.hh>

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

    auto decomposition = Dumux::Mortar::DecompositionFactory<MortarGridGeometry, SubDomainGridGeometry>{}
                            .withMortar(mortarGG)
                            .withSubDomain(topGG)
                            .withSubDomain(bottomGG)
                            .make();

    int exitCode = 0;
    if (decomposition.numberOfMortars() != 1 or decomposition.numberOfSubDomains() != 2)
    {
        std::cout << "Unexpected subdomain/mortar counts" << std::endl;
        exitCode = 1;
    }

    std::vector<int> errors;
    decomposition.visitCoupledMortarsOf(*topGG, [&] (const auto& mortar) { errors.push_back(mortar != mortarGG); });
    decomposition.visitCoupledMortarsOf(*bottomGG, [&] (const auto& mortar) { errors.push_back(mortar != mortarGG); });
    decomposition.visitCoupledSubDomainsOf(*mortarGG, [&] (const auto& subDomain) {
        errors.push_back(subDomain != topGG && subDomain != bottomGG);
    });

    if (std::accumulate(errors.begin(), errors.end(), 0) > 0 || errors.size() != 4)
    {
        std::cout << "Unexpected mappings" << std::endl;
        exitCode = 1;
    }

    errors.clear();
    decomposition.visitSubDomainTraceWith(*mortarGG, *topGG, [&] (const auto& trace) {
        errors.push_back(trace.gridView().size(0) != 10);
    });
    decomposition.visitSubDomainTraceWith(*mortarGG, *bottomGG, [&] (const auto& trace) {
        errors.push_back(trace.gridView().size(0) != 10);
    });
    if (std::accumulate(errors.begin(), errors.end(), 0) > 0 || errors.size() != 2)
    {
        std::cout << "Unexpected traces" << std::endl;
        exitCode = 1;
    }

    return exitCode;
}
