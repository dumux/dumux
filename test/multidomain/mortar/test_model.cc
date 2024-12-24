#include <unordered_map>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/multidomain/mortar/decomposition.hh>
#include <dumux/multidomain/mortar/model.hh>


template<typename MortarSolution, typename MortarGrid, typename GridGeometry>
struct Solver : public Dumux::Mortar::SubDomainSolver<MortarSolution, MortarGrid, GridGeometry>
{
    using ParentType = Dumux::Mortar::SubDomainSolver<MortarSolution, MortarGrid, GridGeometry>;
    using ParentType::ParentType;
    using typename ParentType::Trace;

    void solve() override {
        wasSolved = true;
    }

    virtual void setTraceVariables(std::size_t mortarId, MortarSolution) override {
        traceWasSet[mortarId] = true;
    }

    virtual void registerMortarTrace(std::shared_ptr<const Trace> trace, std::size_t mortarId) override {
        numberOfTraceElements[mortarId] = trace->grid().leafGridView().size(0);
    }

    virtual MortarSolution assembleTraceVariables(std::size_t mortarId) const override {
        MortarSolution result;
        result.resize(numberOfTraceElements.at(mortarId));
        result = 0;
        return result;
    }

    bool wasSolved = false;
    std::unordered_map<std::size_t, bool> traceWasSet;
    std::unordered_map<std::size_t, std::size_t> numberOfTraceElements;
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

    auto model = Dumux::Mortar::ModelFactory<MortarSolution, MortarGridGeometry, SubDomainGridGeometry>{}
        .withMortar(mortarGG)
        .withSubDomain(topSolver)
        .withSubDomain(bottomSolver)
        .make();

    MortarSolution x(model.numMortarDofs());

    int exitCode = 0;
    if (x.size() != mortarGG->numDofs())
    {
        std::cout << "Unexpected mortar solution size" << std::endl;
        exitCode = 1;
    }

    model.solveSubDomains();
    model.setMortar(x);
    model.assembleMortarResidual(x);

    if (not bottomSolver->wasSolved) { std::cout << "Solve was not invoked" << std::endl; exitCode += 1; }
    if (not topSolver->wasSolved) { std::cout << "Solve was not invoked" << std::endl; exitCode += 1; }

    if (bottomSolver->traceWasSet.size() != 1 or not bottomSolver->traceWasSet.at(0))
    { std::cout << "Trace not set" << std::endl; exitCode += 1; }

    if (topSolver->traceWasSet.size() != 1 or not topSolver->traceWasSet.at(0))
    { std::cout << "Trace not set" << std::endl; exitCode += 1; }

    if (bottomSolver->numberOfTraceElements.size() != 1 or bottomSolver->numberOfTraceElements.at(0) != 10)
    { std::cout << "Unexpected trace" << std::endl; exitCode += 1; }
    if (topSolver->numberOfTraceElements.size() != 1 or topSolver->numberOfTraceElements.at(0) != 10)
    { std::cout << "Unexpected trace" << std::endl; exitCode += 1; }

    return exitCode;
}
