#include <unordered_map>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/multidomain/mortar/decomposition.hh>
#include <dumux/multidomain/mortar/interfaceoperator.hh>
#include <dumux/multidomain/mortar/model.hh>


template<typename MortarSolution, typename MortarGrid, typename GridGeometry>
struct Solver : public Dumux::Mortar::SubDomainSolver<MortarSolution, MortarGrid, GridGeometry>
{
    using ParentType = Dumux::Mortar::SubDomainSolver<MortarSolution, MortarGrid, GridGeometry>;
    using ParentType::ParentType;
    using typename ParentType::Trace;

    void solve() override {}
    virtual void setTraceVariables(std::size_t, MortarSolution) override {}

    virtual void registerMortarTrace(std::shared_ptr<const Trace> trace, std::size_t mortarId) override
    { numberOfTraceElements[mortarId] = trace->grid().leafGridView().size(0); }

    virtual MortarSolution assembleTraceVariables(std::size_t mortarId) const override {
        MortarSolution result;
        result.resize(numberOfTraceElements.at(mortarId));
        result = 0;
        return result;
    }

    std::unordered_map<size_t, std::size_t> numberOfTraceElements;
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

    using MortarSolution = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using SubDomainSolver = Solver<MortarSolution, MortarGrid, SubDomainGridGeometry>;
    auto bottomSolver = std::make_shared<SubDomainSolver>(bottomGG);
    auto topSolver = std::make_shared<SubDomainSolver>(topGG);

    Dumux::Mortar::InterfaceOperator op{
        Dumux::Mortar::ModelFactory<MortarSolution, MortarGridGeometry, SubDomainGridGeometry>{}
            .withMortar(mortarGG)
            .withSubDomain(topSolver)
            .withSubDomain(bottomSolver)
            .make()
    };

    MortarSolution x(mortarGG->numDofs());
    MortarSolution r(mortarGG->numDofs());
    x = 0.0;
    op.apply(x, r);
    op.applyscaleadd(0.5, x, r);

    return 0;
}
