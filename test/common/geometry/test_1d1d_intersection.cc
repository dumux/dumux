#include <config.h>

#include <iostream>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/geometry/geometryintersection.hh>

#ifndef DOXYGEN
template<int dimWorld>
Dune::FieldVector<double, dimWorld>
makePosition( double x, double y )
{
    if (dimWorld == 2)
        return {x, y};
    else if (dimWorld == 3)
        return {x, y, 1.0};
    else
        DUNE_THROW(Dune::InvalidStateException, "Invalid dimworld");
}

template<int dimWorld>
Dune::MultiLinearGeometry<double, 1, dimWorld>
makeLine( Dune::FieldVector<double, dimWorld>&& source,
          Dune::FieldVector<double, dimWorld>&& target)
{
    using CornerStorage = std::vector<Dune::FieldVector<double, dimWorld>>;
    return { Dune::GeometryTypes::line,
             CornerStorage({std::move(source), std::move(target)}) };
}

template<int dimworld>
bool testPointIntersection(const Dune::MultiLinearGeometry<double, 1, dimworld>& seg1,
                           const Dune::MultiLinearGeometry<double, 1, dimworld>& seg2,
                           bool foundExpected)
{
    using Segment = Dune::MultiLinearGeometry<double, 1, dimworld>;
    using Policy = Dumux::IntersectionPolicy::PointPolicy<double, dimworld>;
    using Algorithm = Dumux::GeometryIntersection<Segment, Segment, Policy>;

    typename Algorithm::Intersection is;
    bool found = Algorithm::intersection(seg1, seg2, is);

    if (!found && foundExpected)
        std::cerr << "Failed detecting intersection with ";
    else if (found && foundExpected)
        std::cout << "Found intersection at " << is << " with ";
    else if (found && !foundExpected)
        std::cerr << "Found false positive: intersection at " << is << " with ";
    else if (!found && !foundExpected)
        std::cout << "No intersection with ";

    std::cout << seg1.corner(0) << " - " << seg1.corner(1) << " / "
              << seg2.corner(0) << " - " << seg2.corner(1) << std::endl;

    return (found == foundExpected);
}

template<int dimWorld>
void testPointIntersections(std::vector<bool>& returns)
{
    // test intersection points of segments
    for (auto scaling : {1.0, 1e3, 1e12, 1e-12})
    {
        std::cout << "Test point intersections with scaling " << scaling << std::endl;

        const auto seg1 = makeLine( makePosition<dimWorld>(0.0, 0.0), makePosition<dimWorld>(0.0, scaling*1.0) );
        const auto seg2 = makeLine( makePosition<dimWorld>(0.0, 0.0), makePosition<dimWorld>(scaling*1.0, 0.0) );
        const auto seg3 = makeLine( makePosition<dimWorld>(scaling*0.5, 0.0), makePosition<dimWorld>(scaling*1.0, 0.0) );
        const auto seg4 = makeLine( makePosition<dimWorld>(-1.0*scaling*0.5, scaling*0.5), makePosition<dimWorld>(scaling*1.0, 0.0) );
        const auto seg5 = makeLine( makePosition<dimWorld>(0.0, scaling*1.0), makePosition<dimWorld>(scaling*1.0, 0.0) );
        const auto seg6 = makeLine( makePosition<dimWorld>(-1.0*scaling*0.5, scaling*0.5), makePosition<dimWorld>(0.0, scaling*0.5) );
        const auto seg7 = makeLine( makePosition<dimWorld>(-1.0*scaling*0.5, scaling*0.5), makePosition<dimWorld>(-1.0*scaling*0.1, scaling*0.5) );

        // segments with same orientation touching in source or target or not at all
        const auto seg8 = makeLine( makePosition<dimWorld>(0.0, 0.0), makePosition<dimWorld>(0.0, -1.0*scaling*1.0) );
        const auto seg9 = makeLine( makePosition<dimWorld>(0.0, scaling*1.0), makePosition<dimWorld>(0.0, scaling*2.0) );
        const auto seg10 = makeLine( makePosition<dimWorld>(0.0, -1.0*0.01*scaling), makePosition<dimWorld>(0.0, -1.0*scaling*1.0) );
        const auto seg11 = makeLine( makePosition<dimWorld>(0.0, scaling*1.01), makePosition<dimWorld>(0.0, scaling*2.0) );

        returns.push_back(testPointIntersection(seg1, seg2, true));
        returns.push_back(testPointIntersection(seg1, seg3, false));
        returns.push_back(testPointIntersection(seg1, seg4, true));
        returns.push_back(testPointIntersection(seg1, seg5, true));
        returns.push_back(testPointIntersection(seg1, seg6, true));
        returns.push_back(testPointIntersection(seg1, seg7, false));
        returns.push_back(testPointIntersection(seg1, seg8, true));
        returns.push_back(testPointIntersection(seg1, seg9, true));
        returns.push_back(testPointIntersection(seg1, seg10, false));
        returns.push_back(testPointIntersection(seg1, seg11, false));

        std::cout << std::endl;
    }
}

#endif

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // collect returns to determine exit code
    std::vector<bool> returns;

    // test for dimWorld = 2
    testPointIntersections<2>(returns);

    // TODO: implement and test for dimWorld = 3
    // testPointIntersections<3>(returns);

    // TODO: implement and test segment intersections

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
