//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/io/format.hh>
#include <dumux/geometry/volume.hh>

#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/facecentered/diamond/geometryhelper.hh>

namespace Mock {

template<int dim, int dimworld>
struct GridView
{
    template<int i>
    struct Codim
    {
        struct Entity { using Geometry = Dune::MultiLinearGeometry<double, dim, dimworld>; };
    };

    // dummy
    struct Intersection { using Geometry = int; };

    using ctype = double;
    static constexpr int dimension = dim;
    static constexpr int dimensionworld = dimworld;
};

template<int dim, int dimworld>
struct SubGeoType
{
    struct Traits
    {
        using CornerStorage = std::vector<Dune::FieldVector<double, dimworld>>;
        using Geometry = Dune::MultiLinearGeometry<double, dim, dimworld>;
        using LocalIndexType = int;
    };
};

template<int dim, int dimworld>
using ScvType = SubGeoType<dim, dimworld>;

template<int dim, int dimworld>
using ScvfType = SubGeoType<dim-1, dimworld>;

} // end namespace Mock


template<int dim, int dimworld>
auto referenceGeometryCorners(Dune::GeometryType type)
{
    using C = Dune::FieldVector<double, dimworld>;
    std::vector<C> corners;
    const auto refElement = Dune::referenceElement<double, dim>(type);

    // make this work for dim < dimworld
    for (int i = 0; i < refElement.size(dim); ++i)
    {
        C corner(0.0);
        const auto refCorner = refElement.position(i, dim);
        for (int k = 0; k < dim; ++k)
            corner[k] = refCorner[k];
        corners.emplace_back(std::move(corner));
    }

    return corners;
}

template<int dim, int dimworld>
Dune::MultiLinearGeometry<double, dim, dimworld>
makeGeometry(Dune::GeometryType type, const std::vector<Dune::FieldVector<double, dimworld>>& corners)
{ return { type, corners }; }


template<int dim, int dimworld>
int testVolumeImpl(Dune::GeometryType type, bool distorted)
{
    int result = 0;

    std::cout << "\nTesting geometry type " << type << " in " << dimworld << "d ------> ";

    auto corners = referenceGeometryCorners<dim, dimworld>(type);
    if (distorted) // some linear transformation
        for (int i = 1; i < corners.size(); ++i)
            for (int j = 0; j < dim; ++j)
                corners[i][j] += 0.5*(j+1)*(corners[i][j]-corners[0][j]);

    const auto geo = makeGeometry<dim, dimworld>(type, corners);
    const auto vol0 = Dumux::volume(geo);
    const auto vol1 = Dumux::convexPolytopeVolume(geo);
    const auto vol2 = geo.volume();

    std::cout << Dumux::Fmt::format("integrated: {:.14e}; Dumux::volume: {:.14e}; geometry.volume(): {:.14e}", vol0, vol1, vol2) << std::endl;

    if (std::abs(vol0-vol1) > std::max(vol0, vol1)*1e-14)
    {
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>> FAILED" << std::endl;
        result += 1;
    }

    if (type != Dune::GeometryTypes::pyramid)
    {
        std::cout << "\n  Testing box scvs on this geometry " << std::endl;

        using Helper = Dumux::BoxGeometryHelper<
            Mock::GridView<dim, dimworld>,
            dim,
            Mock::ScvType<dim, dimworld>,
            Mock::ScvfType<dim, dimworld>
        >;
        Helper helper(geo);

        for (int i = 0; i < corners.size(); ++i)
        {
            const auto scvCorners = helper.getScvCorners(i);
            const auto geo = makeGeometry<dim, dimworld>(Dune::GeometryTypes::cube(dim), scvCorners);
            const auto vol0 = Dumux::volume(geo);
            const auto vol1 = Dumux::convexPolytopeVolume(geo);
            const auto vol2 = geo.volume();
            std::cout << Dumux::Fmt::format("   - scv {} integrated: {:.14e}; Dumux::volume: {:.14e}; geometry.volume(): {:.14e}", i, vol0, vol1, vol2) << std::endl;
            if (std::abs(vol0-vol1) > std::max(vol0, vol1)*1e-14)
            {
                std::cout << "    >>>>>>>>>>>>>>>>>>>>>>> FAILED" << std::endl;
                result += 1;
            }
        }
    }

    if (dim > 1 && type != Dune::GeometryTypes::pyramid && type != Dune::GeometryTypes::prism)
    {
        std::cout << "\n  Testing diamond scvs on this geometry " << std::endl;

        using Helper = Dumux::DiamondGeometryHelper<
            Mock::GridView<dim, dimworld>,
            Mock::ScvType<dim, dimworld>,
            Mock::ScvfType<dim, dimworld>
        >;
        Helper helper(geo);

        const auto scvType = type == Dune::GeometryTypes::hexahedron ? Dune::GeometryTypes::pyramid : Dune::GeometryTypes::simplex(dim);

        for (int i = 0; i < helper.numScv(); ++i)
        {
            const auto scvCorners = helper.getScvCorners(i);
            const auto geo = makeGeometry<dim, dimworld>(scvType, scvCorners);
            const auto vol0 = Dumux::volume(geo);
            const auto vol1 = Dumux::convexPolytopeVolume(geo);
            const auto vol2 = geo.volume();
            std::cout << Dumux::Fmt::format("   - scv {} integrated: {:.14e}; Dumux::volume: {:.14e}; geometry.volume(): {:.14e}", i, vol0, vol1, vol2) << std::endl;
            if (std::abs(vol0-vol1) > std::max(vol0, vol1)*1e-14)
            {
                std::cout << "    >>>>>>>>>>>>>>>>>>>>>>> FAILED" << std::endl;
                result += 1;
            }
        }
    }

    return result;
}

template<int dim, int dimworld>
int testVolumeDistorted(Dune::GeometryType type)
{ return testVolumeImpl<dim, dimworld>(type, true); }

template<int dim, int dimworld>
int testVolume(Dune::GeometryType type)
{ return testVolumeImpl<dim, dimworld>(type, false); }

int main(int argc, char* argv[])
{
    int result = 0;
    std::cout << "Tests on reference element " << std::endl;
    result += testVolume<1, 1>(Dune::GeometryTypes::simplex(1));
    result += testVolume<1, 2>(Dune::GeometryTypes::simplex(1));
    result += testVolume<1, 3>(Dune::GeometryTypes::simplex(1));
    result += testVolume<2, 2>(Dune::GeometryTypes::simplex(2));
    result += testVolume<2, 3>(Dune::GeometryTypes::simplex(2));
    result += testVolume<3, 3>(Dune::GeometryTypes::simplex(3));
    result += testVolume<2, 2>(Dune::GeometryTypes::cube(2));
    result += testVolume<2, 3>(Dune::GeometryTypes::cube(2));
    result += testVolume<3, 3>(Dune::GeometryTypes::cube(3));
    result += testVolume<3, 3>(Dune::GeometryTypes::pyramid);
    result += testVolume<3, 3>(Dune::GeometryTypes::prism);

    std::cout << "\nTests on distorted element " << std::endl;
    result += testVolumeDistorted<1, 1>(Dune::GeometryTypes::simplex(1));
    result += testVolumeDistorted<1, 2>(Dune::GeometryTypes::simplex(1));
    result += testVolumeDistorted<1, 3>(Dune::GeometryTypes::simplex(1));
    result += testVolumeDistorted<2, 2>(Dune::GeometryTypes::simplex(2));
    result += testVolumeDistorted<2, 3>(Dune::GeometryTypes::simplex(2));
    result += testVolumeDistorted<3, 3>(Dune::GeometryTypes::simplex(3));
    result += testVolumeDistorted<2, 2>(Dune::GeometryTypes::cube(2));
    result += testVolumeDistorted<2, 3>(Dune::GeometryTypes::cube(2));
    result += testVolumeDistorted<3, 3>(Dune::GeometryTypes::cube(3));
    result += testVolumeDistorted<3, 3>(Dune::GeometryTypes::pyramid);
    result += testVolumeDistorted<3, 3>(Dune::GeometryTypes::prism);

    if (result > 0)
    {
        std::cout << "\n###########################################\n"
                  << result << " SUB-TESTS FAILED\n"
                  << "###########################################\n";
    }

    return result;
}
