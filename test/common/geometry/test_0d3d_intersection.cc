#include <config.h>

#include <iostream>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/geometry/intersectspointgeometry.hh>

#ifndef DOXYGEN
template<int dimworld = 3, class Geometry>
bool testIntersection(const Geometry& geo,
                      const Dune::FieldVector<double, dimworld>& p,
                      bool foundExpected = false)
{
    bool found = Dumux::intersectsPointGeometry(p, geo);
    if (!found && foundExpected)
        std::cerr << "Failed detecting intersection with " << p << std::endl;
    else if (found && foundExpected)
        std::cout << "Found intersection with " << p << std::endl;
    else if (found && !foundExpected)
        std::cerr << "Found false positive: intersection with " << p << std::endl;
    else if (!found && !foundExpected)
        std::cout << "No intersection with " << p << std::endl;
    return (found == foundExpected);
}

#endif

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);
    
    constexpr int dimworld = 3;
    
    // collect returns to determine exit code
    std::vector<bool> returns;

    for (auto scaling : {1.0, 1e3, 1e12, 1e-12})
    {
        std::cout << "Test with scaling " << scaling << std::endl;

        using Points = std::vector<Dune::FieldVector<double, dimworld>>;
        using Geo = Dune::MultiLinearGeometry<double, 3, dimworld>;
    
        // test tetrahedron-point intersections
        std::cout << "test tetrahedron-point intersections" << std::endl;

        auto cornersTetrahedron = Points({{0.0, 0.0, 0.0}, {1.0*scaling, 0.0, 0.0}, {0.0, 1.0*scaling, 0.0}, {0.0, 0.0, 1.0*scaling}});
        auto tetrahedron = Geo(Dune::GeometryTypes::tetrahedron, cornersTetrahedron);

        returns.push_back(testIntersection(tetrahedron, {0.0, 0.0, 0.0}, true));
        returns.push_back(testIntersection(tetrahedron, {0.0, 0.0, 0.5*scaling}, true));
        returns.push_back(testIntersection(tetrahedron, {0.25*scaling, 0.25*scaling, 0.5*scaling}, true));
        returns.push_back(testIntersection(tetrahedron, {0.5*scaling, 0.5*scaling, 0.5*scaling}));
        returns.push_back(testIntersection(tetrahedron, {1.01*scaling, 0.0, 0.0}));
        returns.push_back(testIntersection(tetrahedron, {0.5*scaling, 0.0, 0.51*scaling}));

        // test hexahedron-point intersections
        std::cout << "test hexahedron-point intersections" << std::endl;

        auto cornersHexahedron = Points({{0.0, 0.0, 0.0}, {1.0*scaling, 0.0, 0.0}, {0.0, 1.0*scaling, 0.0}, {1.0*scaling, 1.0*scaling, 0.0}, {0.0, 0.0, 1.0*scaling}, {1.0*scaling, 0.0, 1.0*scaling}, {0.0, 1.0*scaling, 1.0*scaling}, {1.0*scaling, 1.0*scaling, 1.0*scaling}});
        auto hexahedron = Geo(Dune::GeometryTypes::hexahedron, cornersHexahedron);

        returns.push_back(testIntersection(hexahedron, {0.0, 0.0, 0.0}, true));

        // test pyramid-point intersections
        std::cout << "test pyramid-point intersections" << std::endl;

        auto cornersPyramid = Points({{0.0, 0.0, 0.0}, {1.0*scaling, 0.0, 0.0}, {0.0, 1.0*scaling, 0.0}, {1.0*scaling, 1.0*scaling, 0.0}, {0.5*scaling, 0.5*scaling, 1.0*scaling}});
        auto pyramid = Geo(Dune::GeometryTypes::pyramid, cornersPyramid);

        returns.push_back(testIntersection(pyramid, {0.0, 0.0, 0.0}, true));

        // test prism-point intersections
        std::cout << "test prism-point intersections" << std::endl;

        auto cornersPrism = Points({{0.0, 0.0, 0.0}, {1.0*scaling, 0.0, 0.0}, {0.0, 1.0*scaling, 0.0}, {0.0, 0.0, 1.0*scaling}, {1.0*scaling, 0.0, 1.0*scaling}, {0.0, 1.0*scaling, 1.0*scaling}});
        auto prism = Geo(Dune::GeometryTypes::prism, cornersPrism);

        returns.push_back(testIntersection(prism, {0.0, 0.0, 0.0}, true));
    }

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](bool i){ return !i; }))
        return 1;

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
