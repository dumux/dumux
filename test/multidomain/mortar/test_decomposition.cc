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

    std::vector<int> exitCodes;
    decomposition.visitCoupledMortarsOf(*topGG, [&] (const auto& mortar) { exitCodes.push_back(mortar != mortarGG); });
    decomposition.visitCoupledMortarsOf(*bottomGG, [&] (const auto& mortar) { exitCodes.push_back(mortar != mortarGG); });
    decomposition.visitCoupledSubDomainsOf(*mortarGG, [&] (const auto& subDomain) {
        exitCodes.push_back(subDomain != topGG && subDomain != bottomGG);
    });

    const auto ec = std::accumulate(exitCodes.begin(), exitCodes.end(), 0);
    if (ec > 0 || exitCodes.size() != 4)
    {
        std::cout << "Mapping test did not succeed" << std::endl;
        return 1;
    }

    return 0;
}
