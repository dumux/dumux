#include <iostream>
#include <memory>
#include <tuple>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/multidomain/mortar/decomposition.hh>


using SubDomainGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
using SubDomainGG = Dumux::CCTpfaFVGridGeometry<typename SubDomainGrid::LeafGridView>;

using MortarGrid = Dune::FoamGrid<1, 2>;
using MortarGG = Dumux::CCTpfaFVGridGeometry<typename MortarGrid::LeafGridView>;
using MortarGridFactory = Dune::StructuredGridFactory<MortarGrid>;

template<typename GGT1, typename GGT2>
void testDecomposition(const Dumux::Mortar::Decomposition<GGT1, GGT2>& decomposition)
{
    if (decomposition.numSubDomains() != 4) DUNE_THROW(Dune::InvalidStateException, "Unexpected number of subdomains");
    if (decomposition.numMortars() != 4) DUNE_THROW(Dune::InvalidStateException, "Unexpected number of mortars");

    for (std::size_t i = 0; i < decomposition.numSubDomains(); ++i)
    {
        int count = 0;
        decomposition.visitInterfacesOfSubdomain(i, [&] (const auto&) { ++count; });
        if (count != 2) DUNE_THROW(Dune::InvalidStateException, "Unexpected number of subdomain->mortar visits");
    }

    for (std::size_t i = 0; i < decomposition.numMortars(); ++i)
    {
        int count = 0;
        decomposition.visitInterfacesOfMortar(i, [&] (const auto&) { ++count; });
        if (count != 2) DUNE_THROW(Dune::InvalidStateException, "Unexpected number of mortar->subdomain visits");
    }
}

auto makeMortarGrid(const Dune::FieldVector<double, 2>& lowerLeft,
                    const Dune::FieldVector<double, 2>& upperRight)
{
    Dune::GridFactory<MortarGrid> factory;
    factory.insertVertex(lowerLeft);
    factory.insertVertex(upperRight);
    factory.insertElement(Dune::GeometryTypes::simplex(1), {0, 1});
    return factory.createGrid();
}

int main() {
    SubDomainGrid sd1{{0.0, 0.0}, {1.0, 1.0}, {10, 10}};
    SubDomainGrid sd2{{1.0, 0.0}, {2.0, 1.0}, {10, 10}};
    SubDomainGrid sd3{{1.0, 1.0}, {2.0, 2.0}, {10, 10}};
    SubDomainGrid sd4{{0.0, 1.0}, {1.0, 2.0}, {10, 10}};

    auto mg1 = makeMortarGrid({0.0, 1.0}, {1.0, 1.0});
    auto mg2 = makeMortarGrid({1.0, 1.0}, {2.0, 1.0});
    auto mg3 = makeMortarGrid({1.0, 0.0}, {1.0, 1.0});
    auto mg4 = makeMortarGrid({1.0, 1.0}, {1.0, 2.0});

    const auto makeSDGG = [] (const auto& grid) { return std::make_shared<SubDomainGG>(grid.leafGridView()); };
    const auto makeMGG = [] (const auto& grid) { return std::make_shared<MortarGG>(grid->leafGridView()); };

    testDecomposition(Dumux::Mortar::Decomposition<SubDomainGG, MortarGG>{
        {makeSDGG(sd1), makeSDGG(sd2), makeSDGG(sd3), makeSDGG(sd4)},
        {makeMGG(mg1), makeMGG(mg2), makeMGG(mg3), makeMGG(mg4)}
    });

    return 0;
}
