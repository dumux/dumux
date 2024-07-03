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
    static_assert(Dumux::Mortar::Decomposition<GGT1, GGT2>::numSubDomains == 4);
    static_assert(Dumux::Mortar::Decomposition<GGT1, GGT2>::numMortars == 4);
}

template<typename... GGs>
struct MultiDomainTraits
{
private:
    template<template<std::size_t> typename T, typename Indices>
    struct TupleImpl;

    template<template<std::size_t> typename T, std::size_t... i>
    struct TupleImpl<T, std::index_sequence<i...>>
    : public std::type_identity<std::tuple<T<i>...>> {};

public:
    template<std::size_t i>
    struct SubDomain {
        using GridGeometry = std::tuple_element_t<i, std::tuple<GGs...>>;
        using PtrType = std::shared_ptr<GridGeometry>;
    };

    static constexpr std::size_t numSubDomains = sizeof...(GGs);

    template<template<std::size_t> typename T>
    using Tuple = typename TupleImpl<T, std::make_index_sequence<numSubDomains>>::type;
};

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

    using SubDomainMDTraits = MultiDomainTraits<SubDomainGG, SubDomainGG, SubDomainGG, SubDomainGG>;
    using MortarMDTraits = MultiDomainTraits<MortarGG, MortarGG, MortarGG, MortarGG>;

    {
        std::cout << "Test from lvalues to MDGridGeometries" << std::endl;
        Dumux::MultiDomainFVGridGeometry<SubDomainMDTraits> subDomainGGs{
            sd1.leafGridView(), sd2.leafGridView(), sd3.leafGridView(), sd4.leafGridView()
        };
        Dumux::MultiDomainFVGridGeometry<MortarMDTraits> mortarGGs{
            mg1->leafGridView(), mg2->leafGridView(), mg3->leafGridView(), mg4->leafGridView()
        };
        testDecomposition(Dumux::Mortar::Decomposition{subDomainGGs, mortarGGs});
    }

    {
        std::cout << "Test from rvalues to MDGridGeometries" << std::endl;
        testDecomposition(Dumux::Mortar::Decomposition{
            Dumux::MultiDomainFVGridGeometry<SubDomainMDTraits>{
                sd1.leafGridView(), sd2.leafGridView(), sd3.leafGridView(), sd4.leafGridView()
            },
            Dumux::MultiDomainFVGridGeometry<MortarMDTraits>{
                mg1->leafGridView(), mg2->leafGridView(), mg3->leafGridView(), mg4->leafGridView()
            }
        });
    }

    return 0;
}
