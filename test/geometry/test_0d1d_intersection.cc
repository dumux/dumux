//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <algorithm>
#include <initializer_list>

#include <dune/common/exceptions.hh>
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

        GlobalPosition delta = b - a;
        delta *= 1.5e-7;
        GlobalPosition p7 = a - delta;
        GlobalPosition p8 = b + delta;
        GlobalPosition p9 = a + delta;
        GlobalPosition p10 = b - delta;

        returns.push_back(testIntersection(a, b, p1, true));
        returns.push_back(testIntersection(a, b, p2, true));
        returns.push_back(testIntersection(a, b, p3, true));
        returns.push_back(testIntersection(a, b, p4));
        returns.push_back(testIntersection(a, b, p5));
        returns.push_back(testIntersection(a, b, p6));

        // test cases where the point is just on or just outside the segment
        returns.push_back(testIntersection(a, b, p7));
        returns.push_back(testIntersection(a, b, p8));
        returns.push_back(testIntersection(a, b, p9, true));
        returns.push_back(testIntersection(a, b, p10, true));

        // test segment that is not axis-parallel
        const GlobalPosition a2(scaling);
        const GlobalPosition b2(scaling*2.0);

        GlobalPosition p11 = a2;
        p11[dimWorld-1] += (b2-a2).two_norm()*1.0e-6;
        returns.push_back(testIntersection(a2, b2, p11));

        // test that triggers bug in 1d2d intersection query
        // which has been fixed with !2274 (found a false positive here)
        GlobalPosition p12(0.0); p12[0] = 0.001969*scaling; p12[1] = 0.0004995*scaling;
        GlobalPosition c(0.0); c[1] = 0.0005*scaling;
        GlobalPosition d = c; d[0] = 0.002*scaling;
        returns.push_back(testIntersection(c, d, p12));
    }
}

#endif

int main (int argc, char *argv[])
{
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
