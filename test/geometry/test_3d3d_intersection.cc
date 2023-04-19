//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// graham convex hull test + triangulation
#include <config.h>

#include <fstream>
#include <iostream>
#include <string>

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

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
