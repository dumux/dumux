// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Test for rotation symmetric grid geometry
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/rotationpolicy.hh>
#include <dumux/discretization/rotationsymmetricscv.hh>
#include <dumux/discretization/rotationsymmetricscvf.hh>
#include <dumux/discretization/rotationsymmetricgridgeometrytraits.hh>

namespace Dumux {

template<class GG>
void runTest(const GG& gg, const double refVolume, const double refSurface)
{
    double volume = 0.0;
    double surface = 0.0;
    for (const auto& element : elements(gg.gridView()))
    {
        auto fvGeometry = localView(gg);
        fvGeometry.bind(element);

        for (const auto& scv : scvs(fvGeometry))
            volume += scv.volume();

        for (const auto& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                surface += scvf.area();
    }

    // compare to reference
    if (!Dune::FloatCmp::eq(refVolume, volume))
        DUNE_THROW(Dune::Exception, "Volume not correct! Reference: " << refVolume << ", computed: " << volume);
    if (!Dune::FloatCmp::eq(refSurface, surface))
        DUNE_THROW(Dune::Exception, "Surface not correct! Reference: " << refSurface << ", computed: " << surface);
}

} // end namespace Dumux

int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // test the disc policy
    {
        using Grid = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
        using GGTraits = RotationSymmetricGridGeometryTraits<CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>, RotationPolicy::disc>;
        using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, /*caching=*/false, GGTraits>;
        using GlobalPosition = typename GridGeometry::SubControlVolume::GlobalPosition;

        // make a grid
        const double innerRadius = 0.1;
        const double outerRadius = 1.0;
        GlobalPosition lower(innerRadius);
        GlobalPosition upper(outerRadius);
        std::array<unsigned int, Grid::dimension> els{{10}};
        std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

        // obtain leaf and make GridGeometry
        auto leafGridView = grid->leafGridView();
        GridGeometry gg(leafGridView);
        gg.update();

        // compute the annulus area and the surface
        const double refVolume = M_PI*(outerRadius*outerRadius - innerRadius*innerRadius);
        const double refSurface = 2.0*M_PI*(innerRadius + outerRadius);
        runTest(gg, refVolume, refSurface);

        std::cout << "Successfully tested disc (annulus) policy." << std::endl;

    } // end disc policy

    // test the ball policy
    {
        using Grid = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
        using GGTraits = RotationSymmetricGridGeometryTraits<CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>, RotationPolicy::ball>;
        using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, /*caching=*/false, GGTraits>;
        using GlobalPosition = typename GridGeometry::SubControlVolume::GlobalPosition;

        // make a grid
        const double innerRadius = 0.1;
        const double outerRadius = 1.0;
        GlobalPosition lower(innerRadius);
        GlobalPosition upper(outerRadius);
        std::array<unsigned int, Grid::dimension> els{{10}};
        std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

        // obtain leaf and make GridGeometry
        auto leafGridView = grid->leafGridView();
        GridGeometry gg(leafGridView);
        gg.update();

        // compute the annulus area and the surface
        const double refVolume = 4.0/3.0*M_PI*(outerRadius*outerRadius*outerRadius - innerRadius*innerRadius*innerRadius);
        const double refSurface = 4.0*M_PI*(innerRadius*innerRadius + outerRadius*outerRadius);
        runTest(gg, refVolume, refSurface);

        std::cout << "Successfully tested ball (shell) policy." << std::endl;

    } // end ball policy

    // test the toroid policy
    {
        using Grid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
        using GGTraits = RotationSymmetricGridGeometryTraits<CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>, RotationPolicy::toroid>;
        using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, /*caching=*/false, GGTraits>;
        using GlobalPosition = typename GridGeometry::SubControlVolume::GlobalPosition;

        // make a grid
        const double innerRadius = 0.1;
        const double outerRadius = 1.0;
        GlobalPosition lower({innerRadius, innerRadius});
        GlobalPosition upper({outerRadius, outerRadius});
        std::array<unsigned int, Grid::dimension> els({20, 20});
        std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

        // obtain leaf and make GridGeometry
        auto leafGridView = grid->leafGridView();
        GridGeometry gg(leafGridView);
        gg.update();

        // compute the annulus area and the surface
        const auto centroidRadius = 0.5*(innerRadius + outerRadius);
        const auto side = (outerRadius - innerRadius);
        const double refVolume = side*side*2.0*M_PI*centroidRadius;
        const double refSurface = 4.0*side*2.0*M_PI*centroidRadius;
        runTest(gg, refVolume, refSurface);

        std::cout << "Successfully tested toroid policy." << std::endl;

    } // end toroid policy

    // test the toroid policy for perfect cylinder
    {
        using Grid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
        using GGTraits = RotationSymmetricGridGeometryTraits<CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>, RotationPolicy::toroid>;
        using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, /*caching=*/false, GGTraits>;
        using GlobalPosition = typename GridGeometry::SubControlVolume::GlobalPosition;

        // make a grid
        const double innerRadius = 0.0;
        const double outerRadius = 1.0;
        GlobalPosition lower({innerRadius, innerRadius});
        GlobalPosition upper({outerRadius, outerRadius});
        std::array<unsigned int, Grid::dimension> els({20, 20});
        std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);

        // obtain leaf and make GridGeometry
        auto leafGridView = grid->leafGridView();
        GridGeometry gg(leafGridView);
        gg.update();

        // compute the annulus area and the surface
        const double refVolume = outerRadius*M_PI*outerRadius*outerRadius;
        const double refSurface = 2.0*M_PI*outerRadius + 2.0*M_PI*outerRadius*outerRadius;
        runTest(gg, refVolume, refSurface);

        std::cout << "Successfully tested toroid policy for perfect cylinder." << std::endl;

    } // end toroid policy

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception &e) {

    std::cout << e << std::endl;
    return 1;
}
