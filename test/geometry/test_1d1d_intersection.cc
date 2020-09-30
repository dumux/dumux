#include <config.h>

#include <iostream>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>

#include <dumux/geometry/geometryintersection.hh>

#ifndef DOXYGEN
template<int dimWorld>
Dune::FieldVector<double, dimWorld>
makePosition( double x, double y, double z = 1.0 )
{
    if (dimWorld == 2)
        return {x, y};
    else if (dimWorld == 3)
        return {x, y, z};
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

template<int dimWorld>
std::string toString(const Dune::FieldVector<double, dimWorld>& p)
{
    std::stringstream stream;
    stream << p;
    return stream.str();
}

template<int dimWorld>
std::string toString(const std::array<Dune::FieldVector<double, dimWorld>, 2>& seg)
{
    std::stringstream stream;
    stream << seg[0] << ", " << seg[1];
    return stream.str();
}

template<int dimworld, class Policy>
bool testSegSegIntersection(const Dune::MultiLinearGeometry<double, 1, dimworld>& seg1,
                            const Dune::MultiLinearGeometry<double, 1, dimworld>& seg2,
                            bool foundExpected)
{
    using Segment = Dune::MultiLinearGeometry<double, 1, dimworld>;
    using Algorithm = Dumux::GeometryIntersection<Segment, Segment, Policy>;

    typename Algorithm::Intersection is;
    bool found = Algorithm::intersection(seg1, seg2, is);

    if (!found && foundExpected)
        std::cerr << "Failed detecting intersection with ";
    else if (found && foundExpected)
        std::cout << "Found intersection with corners " << toString(is) << " with ";
    else if (found && !foundExpected)
        std::cerr << "Found false positive: with corners " << toString(is) << " with ";
    else if (!found && !foundExpected)
        std::cout << "No intersection with ";

    std::cout << seg1.corner(0) << " - " << seg1.corner(1) << " / "
              << seg2.corner(0) << " - " << seg2.corner(1) << std::endl;

    return (found == foundExpected);
}

template<int dimworld>
bool testPointIntersection(const Dune::MultiLinearGeometry<double, 1, dimworld>& seg1,
                           const Dune::MultiLinearGeometry<double, 1, dimworld>& seg2,
                           bool foundExpected)
{
    using PP = Dumux::IntersectionPolicy::PointPolicy<double, dimworld>;
    return testSegSegIntersection<dimworld, PP>(seg1, seg2, foundExpected);
}

template<int dimworld>
bool testSegmentIntersection(const Dune::MultiLinearGeometry<double, 1, dimworld>& seg1,
                           const Dune::MultiLinearGeometry<double, 1, dimworld>& seg2,
                           bool foundExpected)
{
    using SP = Dumux::IntersectionPolicy::SegmentPolicy<double, dimworld>;
    return testSegSegIntersection<dimworld, SP>(seg1, seg2, foundExpected);
}

template<int dimWorld>
void testPointIntersections(std::vector<bool>& returns)
{
    // test intersection points of segments
    for (auto scaling : {1.0, 1e3, 1e12, 1e-12})
    {
        std::cout << "Test point intersections in "<< dimWorld <<"d  with scaling " << scaling << std::endl;

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

template<int dimWorld>
void testSegmentIntersections(std::vector<bool>& returns)
{
    // test intersection points of segments
    for (auto s : {1.0, 1e3, 1e12, 1e-12})
    {
        std::cout << "Test segment intersections in "<< dimWorld <<"d  with scaling " << s << std::endl;

        const auto seg1 = makeLine( makePosition<dimWorld>(0, 1.0*s, 0), makePosition<dimWorld>(0, 0, 0) );
        // reverse
        const auto seg2 = makeLine( makePosition<dimWorld>(0, 0, 0), makePosition<dimWorld>(0, 1.0*s, 0) );
        // touching
        const auto seg3 = makeLine( makePosition<dimWorld>(0, 1.0*s, 0), makePosition<dimWorld>(0, 2.0*s, 0) );
        const auto seg4 = makeLine( makePosition<dimWorld>(0, 0, 0), makePosition<dimWorld>(0, -0.1*s, 0) );
        const auto seg5 = makeLine( makePosition<dimWorld>(0, 1.0*s, 0), makePosition<dimWorld>(1.0*s, 1.0*s, 0) );
        // small
        const auto seg6 = makeLine( makePosition<dimWorld>(0, 0.5*s, 0), makePosition<dimWorld>(0, 0.500001*s, 0) );
        // parallel
        const auto seg7 = makeLine( makePosition<dimWorld>(1e-5*s, 1.0*s, 0), makePosition<dimWorld>(1e-5*s, 0, 0) );
        // intersections
        const auto seg8 = makeLine( makePosition<dimWorld>(0, 0, 0), makePosition<dimWorld>(0, 0.5*s, 0) );
        const auto seg9 = makeLine( makePosition<dimWorld>(0, 1.0*s, 0), makePosition<dimWorld>(0, 0.5*s, 0) );
        const auto seg10 = makeLine( makePosition<dimWorld>(0, -1.0*s, 0), makePosition<dimWorld>(0, 0.5*s, 0) );
        const auto seg11 = makeLine( makePosition<dimWorld>(0, 2.0*s, 0), makePosition<dimWorld>(0, 0.5*s, 0) );
        const auto seg12 = makeLine( makePosition<dimWorld>(0, -2.0*s, 0), makePosition<dimWorld>(0, 2.5*s, 0) );
        // smaller or greater & colinear
        const auto seg13 = makeLine( makePosition<dimWorld>(0, -1.0*s, 0), makePosition<dimWorld>(0, -2.0*s, 0) );
        const auto seg14 = makeLine( makePosition<dimWorld>(0, 2.0*s, 0), makePosition<dimWorld>(0, 2.5*s, 0) );
        // degenerated
        const auto seg15 = makeLine( makePosition<dimWorld>(0, 0, 0), makePosition<dimWorld>(0, 0, 0) );
        const auto seg16 = makeLine( makePosition<dimWorld>(0, 0, 0), makePosition<dimWorld>(0, 1e-20, 0) );

        returns.push_back(testSegmentIntersection(seg1, seg2, true));
        returns.push_back(testSegmentIntersection(seg1, seg3, false));
        returns.push_back(testSegmentIntersection(seg1, seg4, false));
        returns.push_back(testSegmentIntersection(seg1, seg5, false));
        returns.push_back(testSegmentIntersection(seg1, seg6, true));
        returns.push_back(testSegmentIntersection(seg1, seg7, false));
        returns.push_back(testSegmentIntersection(seg1, seg8, true));
        returns.push_back(testSegmentIntersection(seg1, seg9, true));
        returns.push_back(testSegmentIntersection(seg1, seg10, true));
        returns.push_back(testSegmentIntersection(seg1, seg11, true));
        returns.push_back(testSegmentIntersection(seg1, seg12, true));
        returns.push_back(testSegmentIntersection(seg1, seg13, false));
        returns.push_back(testSegmentIntersection(seg1, seg14, false));
        returns.push_back(testSegmentIntersection(seg1, seg15, false));
        returns.push_back(testSegmentIntersection(seg15, seg1, false));
        returns.push_back(testSegmentIntersection(seg1, seg16, false));
        returns.push_back(testSegmentIntersection(seg16, seg1, false));

        std::cout << std::endl;
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
    testPointIntersections<2>(returns);

    // TODO: implement and test for dimWorld = 3
    // testPointIntersections<3>(returns);

    // test segment intersection for dimWorld = 2
    testSegmentIntersections<2>(returns);

    // test segment intersection for dimWorld = 3
    testSegmentIntersections<3>(returns);

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](bool i){ return !i; }))
        return 1;

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
