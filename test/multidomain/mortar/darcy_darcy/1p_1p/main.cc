#include <cmath>
#include <fstream>
#include <memory>

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

#include <dumux/multidomain/mortar/interfacesolver.hh>
#include <dumux/multidomain/mortar/model.hh>
#include <dumux/multidomain/mortar/solvers.hh>

#include "../../dumpoperators.hh"
#include "properties.hh"

#if MORTAR_BOX
using TypeTag = Dumux::Properties::TTag::OnePDarcyMortarBox;
#else
using TypeTag = Dumux::Properties::TTag::OnePDarcyMortarTpfa;
#endif
using Grid = Dumux::GetPropType<TypeTag, Dumux::Properties::Grid>;
using GridGeometry = Dumux::GetPropType<TypeTag, Dumux::Properties::GridGeometry>;
using MortarGrid = Dumux::GetPropType<TypeTag, Dumux::Properties::MortarGrid>;
using MortarSolution = Dumux::GetPropType<TypeTag, Dumux::Properties::MortarSolutionVector>;
using MortarGridGeometry = Dumux::BoxFVGridGeometry<double, typename MortarGrid::LeafGridView>;

using SubDomainSolver = Dumux::Mortar::DefaultSubDomainSolver<TypeTag>;

template<typename Geometry>
auto bboxOf(const Geometry& geo)
{
    using ctype = typename Geometry::ctype;
    using V = typename Geometry::GlobalCoordinate;

    std::array<V, 2> bbox;
    std::ranges::fill(bbox[0], std::numeric_limits<ctype>::max());
    std::ranges::fill(bbox[1], std::numeric_limits<ctype>::min());
    for (int c = 0; c < geo.corners(); ++c)
        for (int d = 0; d < Geometry::coorddimension; ++d)
        {
            bbox[0][d] = std::min(bbox[0][d], geo.corner(c)[d]);
            bbox[1][d] = std::max(bbox[1][d], geo.corner(c)[d]);
        }
    return bbox;
}

template<typename V>
auto makeMortarGrid(const std::array<V, 2> boundingBox, std::array<int, 1> numCells) requires(V::size() == 2)
{
    Dune::GridFactory<MortarGrid> factory;
    auto p = boundingBox[0];
    auto dv = boundingBox[1] - p;
    dv *= 1.0/static_cast<double>(numCells[0]);

    factory.insertVertex(p);
    for (unsigned int i = 0; i < numCells[0]; ++i)
    {
        p += dv;
        factory.insertVertex(p);
        factory.insertElement(Dune::GeometryTypes::line, {i, i + 1});
    }

    return factory.createGrid();
}

template<std::size_t dim>
auto makeGrids(std::array<int, dim> subDomains, unsigned int numCellsPerSubdomain = 10)
{
    Dune::FieldVector<typename Grid::ctype, dim> origin(0);
    Dune::FieldVector<typename Grid::ctype, dim> size(1);
    std::array<int, dim> subDomainCells;
    std::array<int, dim-1> mortarCells;
    std::ranges::fill(subDomainCells, numCellsPerSubdomain);
    std::ranges::fill(mortarCells, std::max(unsigned{1}, static_cast<unsigned int>(numCellsPerSubdomain*0.5)));

    std::vector<std::unique_ptr<Grid>> sdGrids;
    std::vector<std::unique_ptr<MortarGrid>> mortarGrids;

    Grid grid{origin, size, subDomains};
    for (const auto& element : elements(grid.leafGridView()))
    {
        auto box = bboxOf(element.geometry());
        sdGrids.push_back(std::make_unique<Grid>(box[0], box[1], subDomainCells));
        for (const auto& is : intersections(grid.leafGridView(), element))
        {
            if (is.boundary())
                continue;
            if (grid.leafGridView().indexSet().index(is.inside())
                > grid.leafGridView().indexSet().index(is.outside()))
                continue;  // do not visit the same facet twice
            mortarGrids.push_back(makeMortarGrid(bboxOf(is.geometry()), mortarCells));
        }
    }

    return std::make_pair(std::move(sdGrids), std::move(mortarGrids));
}

auto makeSolvers(const std::vector<std::unique_ptr<Grid>>& subDomainGrids)
{
    std::vector<std::shared_ptr<SubDomainSolver>> solvers;
    for (const auto& sdGrid : subDomainGrids)
        solvers.push_back(std::make_shared<SubDomainSolver>(
            std::make_shared<GridGeometry>(sdGrid->leafGridView())
        ));
    return solvers;
}

auto makeModel(const std::vector<std::shared_ptr<SubDomainSolver>>& subDomainSolvers,
               const std::vector<std::unique_ptr<MortarGrid>>& mortarGrids)
{
    Dumux::Mortar::ModelFactory<MortarSolution, MortarGridGeometry, GridGeometry> factory;
    for (const auto& sdSolver : subDomainSolvers)
        factory.insertSubDomain(sdSolver);
    for (const auto& mortarGrid : mortarGrids)
        factory.insertMortar(std::make_shared<MortarGridGeometry>(mortarGrid->leafGridView()));
    auto m = factory.make();
    return std::make_shared<decltype(m)>(std::move(m));
}

int main(int argc, char** argv) {
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init(argc, argv);

    #ifndef SUBDOMAINS_2X1
#define SUBDOMAINS_2X1 PERMEABILITY_JUMP
#endif
    const auto cells = Dumux::getParam<unsigned int>("Grid.CellsPerSubDomain", 10);
#if SUBDOMAINS_2X1
    auto [sdGrids, mortarGrids] = makeGrids<2>({2, 1}, cells);
#else
    auto [sdGrids, mortarGrids] = makeGrids<2>({3, 3}, cells);
#endif
    std::cout << "Cells per subdomain: " << cells << std::endl;
    auto solvers = makeSolvers(sdGrids);
    auto model = makeModel(solvers, mortarGrids);

    std::cout << "Mortar degrees of freedom: " << model->numMortarDofs() << std::endl;

    // in the flux variant the mortar carries the interface flux
    const bool fluxVariant = Dumux::getParam<std::string>("Mortar.Variant", "pressure") == "flux";
    if (fluxVariant)
        model->setCouplingMode(Dumux::Mortar::CouplingMode::natural);

    Dumux::Mortar::InterfaceSolver solver{model};
    if (solver.isConstrained())
        std::cout << "Constraining " << solver.constraints().size() << " floating subdomain(s)" << std::endl;

    MortarSolution x(model->numMortarDofs()); x = 0.0;

    // dense dumps of the homogeneous operator and preconditioner for external eigenvalue analysis
    if (Dumux::getParam<bool>("Mortar.DumpOperators", false))
    {
        MortarSolution delta(x.size());
        model->setHomogeneous(false);
        solver.linearOperator().apply(x, delta);
        delta *= -1.0;
        model->setHomogeneous(true);
        using Dumux::Mortar::Test::dumpMortarOperator; using Dumux::Mortar::Test::dumpMortarVector;
        dumpMortarVector("op_rhs.txt", delta);
        dumpMortarOperator<MortarSolution>("op_S.txt", x.size(), [&] (const auto& e, auto& c) { solver.linearOperator().apply(e, c); });
        if (Dumux::getParam<std::string>("Mortar.Preconditioner", "none") == "interface")
            dumpMortarOperator<MortarSolution>("op_B.txt", x.size(), [&] (const auto& e, auto& c) { solver.preconditioner().apply(c, e); });
        std::cout << "DUMPED operators for " << x.size() << " dofs" << std::endl;
        return 0;
    }

    // A preconditioned conjugate gradient method requires the preconditioner to be symmetric
    // and definite with the same sign as the operator. Both are checked on pseudo-random
    // vectors, since a violation of either is silent in the iteration itself.
    if (Dumux::getParam<bool>("Mortar.CheckPreconditioner", false))
    {
        model->setHomogeneous(true);
        MortarSolution u(x.size()), w(x.size()), Bu(x.size()), Bw(x.size()), Su(x.size());
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            u[i] = std::sin(3.0 + 7.0*static_cast<double>(i));
            w[i] = std::cos(1.0 + 3.0*static_cast<double>(i));
        }
        solver.preconditioner().apply(Bu, u);
        solver.preconditioner().apply(Bw, w);
        solver.linearOperator().apply(u, Su);
        const auto uBw = u*Bw, wBu = w*Bu, uBu = u*Bu, uSu = u*Su;
        std::cout << "Preconditioner symmetry defect: " << std::abs(uBw - wBu)
                  << " on values of size " << std::abs(uBw) << std::endl;
        std::cout << "Signs of (u, Bu) and (u, Su): " << uBu << ", " << uSu << std::endl;
        if (std::abs(uBw - wBu) > 1e-8*std::abs(uBw))
            DUNE_THROW(Dune::MathError, "Preconditioner is not symmetric");
        if (uBu*uSu <= 0.0)
            DUNE_THROW(Dune::MathError, "Preconditioner and operator differ in definiteness sign");
    }

    const auto result = solver.solve(x);
    std::cout << "ITERATIONS " << result.iterations << std::endl;
    if (solver.isConstrained())
        std::cout << "Constraint violation of the solution: " << solver.constraints().violation(x) << std::endl;

    // The exact solution is piecewise linear per subdomain and reproduced by two-point flux
    // approximation on an axis-aligned grid, and it lies in the mortar's own P1 space along a
    // straight interface. The mortar carries a pressure, so it must equal that solution at its
    // vertices. With a permeability jump this also tests the flux continuity the method
    // enforces, since pressure and flux are continuous there only for the right interface value.
    {
        // The axis an interface is normal to is the one it has no extent along.
        const auto normalAxis = [] (const auto& mortarGG) -> int {
            return mortarGG.bBoxMax()[0] - mortarGG.bBoxMin()[0] < 1e-12 ? 0 : 1;
        };

        // The exact flux is stated with respect to the positive axis, while the mortar's
        // reference normal points out of the subdomain with orientation sign +1, so the two
        // agree up to the side that subdomain lies on.
        const auto orientationFactor = [&] (const auto& mortarGG) -> double {
            const auto mortarId = model->decomposition().id(mortarGG);
            const auto axis = normalAxis(mortarGG);
            for (const auto& solverPtr : solvers)
                if (model->orientation(model->decomposition().id(*solverPtr->gridGeometry()), mortarId) == 1)
                {
                    const auto& gg = *solverPtr->gridGeometry();
                    const double subDomainCentre = 0.5*(gg.bBoxMin()[axis] + gg.bBoxMax()[axis]);
                    return subDomainCentre < mortarGG.bBoxMin()[axis] ? 1.0 : -1.0;
                }
            DUNE_THROW(Dune::InvalidStateException, "No positively oriented subdomain for this mortar");
        };

        double maxError = 0.0, largest = 0.0;
        model->decomposition().visitMortars([&] (const auto& mortarPtr) {
            const auto p = model->extractEntriesFor(*mortarPtr, x);
            const int axis = normalAxis(*mortarPtr);
            const double factor = fluxVariant ? orientationFactor(*mortarPtr) : 1.0;
            for (const auto& v : vertices(mortarPtr->gridView()))
            {
                const auto pos = v.geometry().center();
                const auto expected = fluxVariant ? factor*Dumux::mortarExactFlux(pos, axis)
                                                  : Dumux::mortarExactSolution(pos);
                const auto value = p[mortarPtr->vertexMapper().index(v)][0];
                maxError = std::max(maxError, std::abs(value - expected));
                largest = std::max(largest, std::abs(expected));
            }
        });
        std::cout << "ERROR " << maxError << std::endl;
        std::cout << "Mortar solution against the exact solution: max error = " << maxError
                  << " on values up to " << largest << std::endl;
        if (largest < 1e-12)
            DUNE_THROW(Dune::InvalidStateException, "The exact solution vanishes on every mortar");
        if (maxError > Dumux::getParam<double>("Mortar.MaxError", 1e-7))
            DUNE_THROW(Dune::MathError, "Mortar solution is off by " << maxError);
    }

    std::cout << "Writing VTK output" << std::endl;
    unsigned int i = 0; for (const auto& solverPtr : solvers)
    {
        Dune::VTKWriter writer{solverPtr->problem().gridGeometry().gridView()};
        if constexpr (Dumux::DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>)
            writer.addVertexData(solverPtr->solution(), "x");
        else
            writer.addCellData(solverPtr->solution(), "x");
        writer.write("subdomain_" + std::to_string(i++));
    }
    i = 0; model->decomposition().visitMortars([&] (const auto& mortarPtr) {
        Dune::VTKWriter writer{mortarPtr->gridView()};
        auto p = model->extractEntriesFor(*mortarPtr, x);
        writer.addVertexData(p, "p");
        writer.write("mortar_" + std::to_string(i++));
    });

    return 0;
}
