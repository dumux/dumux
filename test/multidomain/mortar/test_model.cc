#include <iostream>
#include <memory>
#include <tuple>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/mortar/model.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "traceproblem.hh"
#include "tracespatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct Mortar { using InheritsFrom = std::tuple<OneP>; };
struct MortarModelTpfa { using InheritsFrom = std::tuple<Mortar, CCTpfaModel>; };
}  // namespace TTAG

template<class TypeTag>
struct Problem<TypeTag, TTag::Mortar>
{ using type = TraceTestProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Mortar>
{
    using type = TraceTestSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Scalar>
    >;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Mortar>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> > ;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::Mortar>
{ using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>; };

}  // namespace Dumux::Properties

using TypeTag = Dumux::Properties::TTag::TYPETAG;

using SubDomainGrid = Dumux::GetPropType<TypeTag, Dumux::Properties::Grid>;
using SubDomainGridGeometry = Dumux::GetPropType<TypeTag, Dumux::Properties::GridGeometry>;
using SubDomainGridVariables = Dumux::GetPropType<TypeTag, Dumux::Properties::GridVariables>;
using SubDomainProblem = Dumux::GetPropType<TypeTag, Dumux::Properties::Problem>;

using MortarGrid = Dune::FoamGrid<1, 2>;
using MortarGridGeometry = Dumux::CCTpfaFVGridGeometry<typename MortarGrid::LeafGridView>;
using MortarGridFactory = Dune::StructuredGridFactory<MortarGrid>;

using SolutionVector = Dune::BlockVector<Dune::FieldVector<double, 1>>;

auto makeMortarGrid(const Dune::FieldVector<double, 2>& lowerLeft,
                    const Dune::FieldVector<double, 2>& upperRight)
{
    Dune::GridFactory<MortarGrid> factory;
    factory.insertVertex(lowerLeft);
    factory.insertVertex(upperRight);
    factory.insertElement(Dune::GeometryTypes::simplex(1), {0, 1});
    return factory.createGrid();
}

class SubDomain
{
 public:
    SubDomain(std::shared_ptr<SubDomainGridVariables> gv)
    : gv_{gv}
    {}

    void solve()
    { DUNE_THROW(Dune::NotImplemented, "solve"); }

    void setMortarTrace(std::size_t, SolutionVector&&)
    { DUNE_THROW(Dune::NotImplemented, "setMortarTrace"); }

    SolutionVector assembleTraceWith(std::size_t i) const
    { DUNE_THROW(Dune::NotImplemented, "assembleTraceWith"); }

    const SubDomainGridGeometry& gridGeometry() const
    { return gv_->gridGeometry(); }

    const std::shared_ptr<SubDomainGridVariables>& gridVariables() const
    { return gv_; }

 private:
    std::shared_ptr<SubDomainGridVariables> gv_;
};

int main(int argc, char** argv) {
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init(argc, argv);

    SubDomainGrid sd1{{0.0, 0.0}, {1.0, 1.0}, {10, 10}};
    SubDomainGrid sd2{{1.0, 0.0}, {2.0, 1.0}, {10, 10}};
    auto mg = makeMortarGrid({1.0, 0.0}, {1.0, 1.0});

    auto sd1GG = std::make_shared<SubDomainGridGeometry>(sd1.leafGridView());
    auto sd2GG = std::make_shared<SubDomainGridGeometry>(sd2.leafGridView());
    auto mGG = std::make_shared<MortarGridGeometry>(mg->leafGridView());

    auto sd1GV = std::make_shared<SubDomainGridVariables>(std::make_shared<SubDomainProblem>(sd1GG), sd1GG);
    auto sd2GV = std::make_shared<SubDomainGridVariables>(std::make_shared<SubDomainProblem>(sd2GG), sd2GG);

    Dumux::Mortar::ModelFactory<SolutionVector, MortarGridGeometry> factory;
    factory.insertMortar(mGG);
    factory.insertSubDomain(std::make_shared<SubDomain>(sd1GV), [&] (auto& subDomain, auto&& trace, std::size_t mortarId) {
        auto op = Dumux::TraceOperator{std::move(trace), [] (auto&&...) { return 1.0; }};
    });

    return 0;
}
