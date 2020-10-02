#include <config.h>

#include <iostream>
#include <algorithm>
#include <initializer_list>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dumux/geometry/intersectspointgeometry.hh>

#ifndef DOXYGEN
template<int dimworld = 3>
bool testIntersection(const Dune::FieldVector<double, dimworld>& a,
                      const Dune::FieldVector<double, dimworld>& b,
                      const Dune::FieldVector<double, dimworld>& p,
                      bool foundExpected = false)
{
    bool found = Dumux::intersectsPointSimplex(p, a, b);
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

template<int dimWorld>
void testIntersections(std::vector<bool>& returns)
{
    // test if points lie on 3d segments
    using GlobalPosition = Dune::FieldVector<double, dimWorld>;

    for (auto scaling : {1.0, 1e3, 1e12, 1e-12})
    {
        const GlobalPosition a(0.0);
        auto b = a;
        b[dimWorld-1] = 1.0*scaling;

        GlobalPosition p1 = a;
        GlobalPosition p2 = b;
        GlobalPosition p3 = a + b;
        p3 /= 2.0;
        GlobalPosition p4 = b + (b - a);
        GlobalPosition p5 = a - (b - a);
        GlobalPosition p6(0.5*scaling);

        returns.push_back(testIntersection(a, b, p1, true));
        returns.push_back(testIntersection(a, b, p2, true));
        returns.push_back(testIntersection(a, b, p3, true));
        returns.push_back(testIntersection(a, b, p4));
        returns.push_back(testIntersection(a, b, p5));
        returns.push_back(testIntersection(a, b, p6));
    }
}

#endif

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // collect returns to determine exit code
    std::vector<bool> returns;

    // test for dimWorld = 2
    testIntersections<2>(returns);

    // test for dimWorld = 3
    testIntersections<3>(returns);

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](bool i){ return !i; }))
        return 1;

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
