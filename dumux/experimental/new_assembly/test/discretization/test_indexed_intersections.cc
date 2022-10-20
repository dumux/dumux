#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <numbers>

#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/indexedintersections.hh>

struct ReferenceFace
{
    Dune::FieldVector<double, 2> center;
    Dune::FieldVector<double, 2> normal;

    template<typename GridImp, typename ISImp>
    bool operator==(const Dune::Intersection<GridImp, ISImp>& is) const
    {
        return Dune::FloatCmp::eq(center, is.geometry().center(), 1e-6)
            && Dune::FloatCmp::eq(normal, is.centerUnitOuterNormal(), 1e-6);
    }
};


template<typename GridImp, typename ISImp>
void removeMatch(const Dune::Intersection<GridImp, ISImp>& is, std::vector<ReferenceFace>& expectedFaces)
{
    auto it = std::ranges::find_if(expectedFaces, [&] (const ReferenceFace& face) { return face == is; });
    if (it == std::ranges::end(expectedFaces))
        DUNE_THROW(Dune::InvalidStateException,
                    "Unexpected face: "
                        << "c = " << is.geometry().center() << ", "
                        << "n = " << is.centerUnitOuterNormal());
    expectedFaces.erase(it);
}

void makeUnique(std::vector<int>& values)
{
    std::ranges::sort(values);
    values.erase(std::unique(values.begin(), values.end()), values.end());
}


void testConforming()
{
    std::cout << "Running conforming test..." << std::flush;

    std::vector<ReferenceFace> expectedFaces{
        {.center = {0.0, 0.5}, .normal = {-1.0, 0.0}},
        {.center = {0.5, 0.0}, .normal = {0.0, -1.0}},
        {.center = {1.0, 0.5}, .normal = {1.0, 0.0}},
        {.center = {0.5, 1.0}, .normal = {0.0, 1.0}}
    };

    std::vector<int> indices;
    Dune::YaspGrid<2> grid{{1.0, 1.0}, {1, 1}};
    for (const auto& element : elements(grid.leafGridView()))
        for (const auto& [i, is] : Dumux::indexedIntersections(grid.leafGridView(), element))
        {
            indices.push_back(i);
            removeMatch(is, expectedFaces);
        }

    if (expectedFaces.size() != 0)
        DUNE_THROW(Dune::InvalidStateException, "Did not visit all expected faces");
    if (indices.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of face indices");
    makeUnique(indices);
    if (indices.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "Non-unique face indices");

    std::cout << " passed" << std::endl;
}

void testNonConforming()
{
#if HAVE_DUNE_ALUGRID
    std::cout << "Running nonconforming test..." << std::flush;

    std::vector<ReferenceFace> expectedFaces{
        {.center = {0.0, 0.5}, .normal = {-1.0, 0.0}},
        {.center = {0.5, 0.0}, .normal = {0.0, -1.0}},
        {.center = {1.0, 0.25}, .normal = {1.0, 0.0}},
        {.center = {1.0, 0.75}, .normal = {1.0, 0.0}},
        {.center = {0.5, 1.0}, .normal = {0.0, 1.0}},

        {.center = {1.0, 0.25}, .normal = {-1.0, 0.0}},
        {.center = {1.25, 0.0}, .normal = {0.0, -1.0}},
        {.center = {1.5, 0.25}, .normal = {1.0, 0.0}},
        {.center = {1.25, 0.5}, .normal = {0.0, 1.0}},

        {.center = {1.5, 0.25}, .normal = {-1.0, 0.0}},
        {.center = {1.75, 0.0}, .normal = {0.0, -1.0}},
        {.center = {2.0, 0.25}, .normal = {1.0, 0.0}},
        {.center = {1.75, 0.5}, .normal = {0.0, 1.0}},

        {.center = {1.0, 0.75}, .normal = {-1.0, 0.0}},
        {.center = {1.25, 0.5}, .normal = {0.0, -1.0}},
        {.center = {1.5, 0.75}, .normal = {1.0, 0.0}},
        {.center = {1.25, 1.0}, .normal = {0.0, 1.0}},

        {.center = {1.5, 0.75}, .normal = {-1.0, 0.0}},
        {.center = {1.75, 0.5}, .normal = {0.0, -1.0}},
        {.center = {2.0, 0.75}, .normal = {1.0, 0.0}},
        {.center = {1.75, 1.0}, .normal = {0.0, 1.0}}
    };

    using Grid = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>;
    Dune::StructuredGridFactory<Grid> factory;
    auto grid = factory.createCubeGrid(
        {0.0, 0.0},
        {2.0, 1.0},
        {2, 1}
    );

    grid->preAdapt();
    for (const auto& e : elements(grid->leafGridView()))
        if (e.geometry().center()[0] > 1.0)
            grid->mark(1, e);
    grid->adapt();
    grid->postAdapt();

    for (const auto& element : elements(grid->leafGridView()))
    {
        std::vector<int> indices;
        for (const auto& [i, is] : Dumux::indexedIntersections(grid->leafGridView(), element))
        {
            indices.push_back(i);
            removeMatch(is, expectedFaces);
        }

        const auto indicesSize = indices.size();
        if (indicesSize != 4 && indicesSize != 5)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected number of face indices");
        makeUnique(indices);
        if (indices.size() != indicesSize)
            DUNE_THROW(Dune::InvalidStateException, "Non-unique face indices");
    }

    std::cout << " passed" << std::endl;
#else
    std::cout << "Skipping non-conforming test because ALUGrid was not found" << std::endl;
#endif
}

void testNetwork()
{
#if HAVE_DUNE_FOAMGRID
    std::cout << "Running test on network grid..." << std::flush;

    using std::numbers::sqrt2;
    std::vector<ReferenceFace> expectedFaces{
        {.center = {-0.5, -0.5}, .normal = {-1.0/sqrt2, -1.0/sqrt2}},
        {.center = {0.0, 0.0}, .normal = {1.0/sqrt2, 1.0/sqrt2}},
        {.center = {0.0, 0.0}, .normal = {1.0/sqrt2, 1.0/sqrt2}},

        {.center = {0.0, 0.0}, .normal = {-1.0/sqrt2, 1.0/sqrt2}},
        {.center = {0.0, 0.0}, .normal = {-1.0/sqrt2, 1.0/sqrt2}},
        {.center = {0.5, -0.5}, .normal = {1.0/sqrt2, -1.0/sqrt2}},

        {.center = {0.0, 0.5}, .normal = {0.0, 1.0}},
        {.center = {0.0, 0.0}, .normal = {0.0, -1.0}},
        {.center = {0.0, 0.0}, .normal = {0.0, -1.0}}
    };

    Dune::GridFactory<Dune::FoamGrid<1, 2>> factory;
    factory.insertVertex({-0.5, -0.5});
    factory.insertVertex({0.0, 0.0});
    factory.insertVertex({0.5, -0.5});
    factory.insertVertex({0.0, 0.5});
    factory.insertElement(Dune::GeometryTypes::line, {0, 1});
    factory.insertElement(Dune::GeometryTypes::line, {1, 2});
    factory.insertElement(Dune::GeometryTypes::line, {1, 3});

    auto grid = factory.createGrid();
    for (const auto& element : elements(grid->leafGridView()))
    {
        std::vector<int> indices;
        for (const auto& [i, is] : Dumux::indexedIntersections(grid->leafGridView(), element))
        {
            indices.push_back(i);
            removeMatch(is, expectedFaces);
        }

        if (indices.size() != 3)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected number of face indices");
        makeUnique(indices);
        if (indices.size() != 2)
            DUNE_THROW(Dune::InvalidStateException, "Non-unique face indices");
    }

    std::cout << " passed" << std::endl;
#else
    std::cout << "Skipping network grid test because FoamGrid was not found" << std::endl;
#endif
}


int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    testConforming();
    testNonConforming();
    testNetwork();

    std::cout << "All tests passed" << std::endl;
    return 0;
}
