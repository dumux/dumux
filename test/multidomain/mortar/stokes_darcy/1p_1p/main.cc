#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/multidomain/mortar/interfaceoperator.hh>
#include <dumux/multidomain/mortar/model.hh>
#include <dumux/multidomain/mortar/solvers.hh>
#include <dumux/multidomain/mortar/preconditioners.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/porousmediumflow/velocityoutput.hh>

#include "properties.hh"

int main(int argc, char** argv) {
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using DarcyTypeTag = Properties::TTag::OnePDarcyMortarTpfa;
    using StokesMomentumTypeTag = Properties::TTag::OnePStokesMortarMomentum;
    using StokesMassTypeTag = Properties::TTag::OnePStokesMortarMass;

    using DarcyGrid = GetPropType<DarcyTypeTag, Properties::Grid>;
    using StokesGrid = GetPropType<StokesMassTypeTag, Properties::Grid>;
    using MortarGrid = GetPropType<DarcyTypeTag, Properties::MortarGrid>;
    DarcyGrid darcyGrid{{0.0, -1.0}, {1.0, 0.0}, {10, 10}};
    StokesGrid stokesGrid{{0.0, 0.0}, {1.0, 1.0}, {10, 10}};
    auto mortarGrid = Dune::StructuredGridFactory<MortarGrid>::createCubeGrid({0.0, 0.0}, {1.0, 0.0}, {5});

    using DarcyGridGeometry = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    using MassGridGeometry = GetPropType<StokesMassTypeTag, Properties::GridGeometry>;
    using MomentumGridGeometry = GetPropType<StokesMomentumTypeTag, Properties::GridGeometry>;
    auto darcyGridGeometry = std::make_shared<DarcyGridGeometry>(darcyGrid.leafGridView());
    auto stokesMassGridGeometry = std::make_shared<MassGridGeometry>(stokesGrid.leafGridView());
    auto stokesMomentumGridGeometry = std::make_shared<MomentumGridGeometry>(stokesGrid.leafGridView());

    using DarcySolver = Mortar::DefaultSubDomainSolver<DarcyTypeTag>;
    using StokesSolver = Mortar::DefaultStaggeredFreeFlowSubDomainSolver<StokesMomentumTypeTag, StokesMassTypeTag>;
    auto darcySolver = std::make_shared<DarcySolver>(darcyGridGeometry);
    auto stokesSolver = std::make_shared<StokesSolver>(stokesMomentumGridGeometry, stokesMassGridGeometry);

    using MortarGridGeometry = BoxFVGridGeometry<double, typename MortarGrid::LeafGridView>;
    auto mortarGridGeometry = std::make_shared<MortarGridGeometry>(mortarGrid->leafGridView());

    using MortarSolution = GetPropType<DarcyTypeTag, Properties::MortarSolutionVector>;
    using ModelFactory = Dumux::Mortar::ModelFactory<MortarSolution, MortarGridGeometry, DarcyGridGeometry, MomentumGridGeometry>;
    auto model = [&] () {
        auto m = ModelFactory{}
            .withSubDomain(stokesSolver)
            .withSubDomain(darcySolver)
            .withMortar(mortarGridGeometry)
            .make();
        return std::make_shared<decltype(m)>(std::move(m));
    } ();

    MortarSolution x(model->numMortarDofs()); x = 0.0;
    auto delta = x;

    std::cout << "Compute initial jump" << std::endl;
    darcySolver->problem().useHomogeneousBoundaryCondition(false);
    stokesSolver->momentumProblem().useHomogeneousBoundaryCondition(false);
    Dumux::Mortar::InterfaceOperator op{model};
    op.apply(x, delta);

    std::cout << "Solve homogeneous problem" << std::endl;
    darcySolver->problem().useHomogeneousBoundaryCondition(true);
    stokesSolver->momentumProblem().useHomogeneousBoundaryCondition(true);
    delta *= -1.0;
    Dune::InverseOperatorResult result;
    Dumux::Mortar::NoPreconditioner<MortarSolution> prec;
    Dune::RestartedGMResSolver<MortarSolution>(op, prec, 1e-8, 5, 1000, 2).apply(x, delta, result);
    if (!result.converged)
        DUNE_THROW(Dune::InvalidStateException, "Linear solver did not converge");

    std::cout << "Compute solution" << std::endl;
    darcySolver->problem().useHomogeneousBoundaryCondition(false);
    stokesSolver->momentumProblem().useHomogeneousBoundaryCondition(false);
    op.apply(x, delta);


    std::cout << "Writing VTK output" << std::endl;
    VtkOutputModule stokesVTKWriter{
        stokesSolver->massGridVariables(),
        stokesSolver->massSolution(),
        "stokes"
    };
    GetPropType<StokesMassTypeTag, Properties::IOFields>::initOutputModule(stokesVTKWriter);
    stokesVTKWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<typename StokesSolver::MassGridVariables>>());
    stokesVTKWriter.write(0.0);

    VtkOutputModule darcyVTKWriter{
        darcySolver->gridVariables(),
        darcySolver->solution(),
        "darcy"
    };
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVTKWriter);
    darcyVTKWriter.addVelocityOutput(
        std::make_shared<
            PorousMediumFlowVelocityOutput<
                typename DarcySolver::GridVariables,
                GetPropType<DarcyTypeTag, Properties::FluxVariables>
            >
        >(darcySolver->gridVariables())
    );
    darcyVTKWriter.write(0.0);

    Dune::VTKWriter writer{mortarGridGeometry->gridView()};
    writer.addVertexData(x, "p");
    writer.write("mortar");

    return 0;
}
