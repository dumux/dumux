//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// graham convex hull test + triangulation
#include <config.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/triangulation.hh>

#include "transformation.hh"
#include "writetriangulation.hh"

template<typename Transformation>
auto makeCube(const Transformation& transform)
{
    using CubeGeometry = Dune::MultiLinearGeometry<double, 3, 3>;
    using GlobalPosition = typename CubeGeometry::GlobalCoordinate;

    return CubeGeometry(
        Dune::GeometryTypes::cube(3),
        std::vector<GlobalPosition>({
            transform({0.0, 0.0, 0.0}), transform({1.0, 0.0, 0.0}),
            transform({0.0, 1.0, 0.0}), transform({1.0, 1.0, 0.0}),
            transform({0.0, 0.0, 1.0}), transform({1.0, 0.0, 1.0}),
            transform({0.0, 1.0, 1.0}), transform({1.0, 1.0, 1.0})
        })
    );
}

// The intersections of a tetrahedron with the cells of a grid covering it partition the
// tetrahedron, so their volumes add up to the tetrahedron volume, also far from the origin.
void testTetrahedronGridPartition()
{
    using Vec = Dune::FieldVector<double, 3>;
    using Geometry = Dune::MultiLinearGeometry<double, 3, 3>;
    using Test = Dumux::GeometryIntersection<Geometry, Geometry>;

    const auto volume = [](const auto& a, const auto& b, const auto& c, const auto& d)
    { return std::abs(Dumux::crossProduct(b - a, c - a)*(d - a))/6.0; };

    const int cells = 24;
    const double h = 1.0/cells;
    for (const double shift : {0.0, 1e3, 1e5})
    {
        std::vector<Vec> tet({
            {0.24019237886466849, 0.67320508075688779, 0.49999999999999994},
            {0.35107743534863856, 0.70928672957597105, 0.54502908912073522},
            {0.30707015948691208, 0.77223143156633334, 0.5366898392529551},
            {0.34212426246196426, 0.74118984276951527, 0.47524190555810586}
        });
        for (auto& corner : tet)
            corner += Vec(shift);
        const Geometry tetrahedron(Dune::GeometryTypes::tetrahedron, tet);

        double sum = 0.0;
        for (int k = 0; k < cells; ++k)
            for (int j = 0; j < cells; ++j)
                for (int i = 0; i < cells; ++i)
                {
                    std::vector<Vec> hexCorners;
                    for (int c = 0; c < 8; ++c)
                        hexCorners.push_back(Vec({shift + (i + (c&1))*h, shift + (j + ((c>>1)&1))*h, shift + (k + ((c>>2)&1))*h}));
                    const Geometry hexahedron(Dune::GeometryTypes::hexahedron, hexCorners);

                    typename Test::Intersection intersection;
                    if (Test::intersection(tetrahedron, hexahedron, intersection))
                        for (const auto& t : Dumux::triangulate<3, 3>(intersection))
                            sum += volume(t[0], t[1], t[2], t[3]);
                }

        const double expected = volume(tet[0], tet[1], tet[2], tet[3]);
        std::cout << "Tetrahedron pieces add up to " << sum << ", tetrahedron volume " << expected << " (shift " << shift << ")" << std::endl;
        if (std::abs(sum - expected) > 1e-10*expected)
            DUNE_THROW(Dune::InvalidStateException, "Volumes of the intersections of a tetrahedron with a grid add up to "
                                                     << sum << " instead of " << expected << " (shift " << shift << ")");
    }
}

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using Vec = Dune::FieldVector<double, 3>;

    int count = 0;
    for (auto scale : {1.0, 1e3, 1e12, 1e-12})
    {
        std::cout << "Testing scale " << scale << std::endl;

        // unit cube with scale
        const auto cube = makeCube(make3DTransformation<double>(scale, Vec(0.0), Vec(1.0), 0.0, false));
        using Geometry3D = decltype(cube);
        using Test = Dumux::GeometryIntersection<Geometry3D, Geometry3D>;
        typename Test::Intersection intersection;

        // cube-cube intersection
        if (Test::intersection(cube, cube, intersection))
        {
            const auto triangulation = Dumux::triangulate<3, 3>(intersection);
            Dumux::writeVTUTetrahedron(triangulation, "cube_intersections");
            if (triangulation.size() != 24)
                DUNE_THROW(Dune::InvalidStateException, "Found " << triangulation.size() << " instead of 24 intersections!");
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "No intersections found!");

        // transform the second cube
        for (const double translation : {0.0, 0.2})
            for (const double angle : {0.0, 0.2*M_PI, 0.5*M_PI, 0.567576567*M_PI, M_PI})
                for (const auto& rotAxis : {Vec({0.0, 0.0, 1.0}), Vec(std::sqrt(3.0)/3.0), Vec({std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0})})
                    if (++count; Test::intersection(cube, makeCube(make3DTransformation<double>(scale, Vec(translation), rotAxis, angle, true)), intersection))
                    {
                        std::cout << "Found intersection (point cloud with " << intersection.size() << " points)" << std::endl;
                        const auto triangulation = Dumux::triangulate<3, 3>(intersection);
                        const std::string fileName = "cube_intersections-" + std::to_string(count);
                        Dumux::writeVTUTetrahedron(triangulation, fileName);
                        std::cout << "Triangulation has " << triangulation.size() << " elements. Dumped to file " << fileName << ".vtu" << std::endl;
                        // TODO: Find better test to check if this is the correct result
                        if (triangulation.size() > 48)
                            DUNE_THROW(Dune::Exception, "Large number of intersections. Something went wrong?");
                    }
    }

    testTetrahedronGridPartition();

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
