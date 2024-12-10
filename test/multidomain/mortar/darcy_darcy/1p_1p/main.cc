#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

#include <dumux/multidomain/mortar/interfaceoperator.hh>
#include <dumux/multidomain/mortar/model.hh>
#include <dumux/multidomain/mortar/solvers.hh>
#include <dumux/multidomain/mortar/preconditioners.hh>

#include "properties.hh"

using TypeTag = Dumux::Properties::TTag::OnePDarcyMortarTpfa;
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

    auto [sdGrids, mortarGrids] = makeGrids<2>({3, 3});
    auto solvers = makeSolvers(sdGrids);
    auto model = makeModel(solvers, mortarGrids);

    MortarSolution x(model->numMortarDofs()); x = 0.0;
    auto delta = x;

    std::cout << "Compute initial jump" << std::endl;
    for (auto solverPtr : solvers)
        solverPtr->problem().useHomogeneousBoundaryCondition(false);
    Dumux::Mortar::InterfaceOperator op{model};
    op.apply(x, delta);

    std::cout << "Solve homogeneous problem" << std::endl;
    for (auto solverPtr : solvers)
        solverPtr->problem().useHomogeneousBoundaryCondition(true);
    delta *= -1.0;
    Dune::InverseOperatorResult result;
    Dumux::Mortar::NoPreconditioner<MortarSolution> prec;
    Dune::CGSolver<MortarSolution>(op, prec, 1e-8, 1000, 1).apply(x, delta, result);
    if (!result.converged)
        DUNE_THROW(Dune::InvalidStateException, "Linear solver did not converge");

    std::cout << "Compute solution" << std::endl;
    for (auto solverPtr : solvers)
        solverPtr->problem().useHomogeneousBoundaryCondition(false);
    op.apply(x, delta);

    std::cout << "Writing VTK output" << std::endl;
    unsigned int i = 0; for (const auto& solverPtr : solvers)
    {
        Dune::VTKWriter writer{solverPtr->problem().gridGeometry().gridView()};
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
