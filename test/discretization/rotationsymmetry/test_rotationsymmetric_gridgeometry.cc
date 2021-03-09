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
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/extrusion.hh>

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
            volume += GG::Extrusion::volume(scv);

        for (const auto& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                surface += GG::Extrusion::area(scvf);

        // compare volume and integrated volume of one scv
        const auto& scv = *(scvs(fvGeometry).begin());
        const auto scvGeometry = scv.geometry();

        double volScv = 0.0;
        const auto ruleScv = Dune::QuadratureRules<double, GG::GridView::dimension>::rule(scvGeometry.type(), 3);
        for (const auto& qp : ruleScv)
            volScv += qp.weight()*GG::Extrusion::integrationElement(scvGeometry, qp.position());

        if (!Dune::FloatCmp::eq(volScv, GG::Extrusion::volume(scv)))
            DUNE_THROW(Dune::Exception, "Integration not correct! Integrated: " << volScv << ", direct: " << GG::Extrusion::volume(scv));

        // compare area and integration area of one scvf
        const auto& scvf = *(scvfs(fvGeometry).begin());
        const auto scvfGeometry = scvf.geometry();

        double volScvf = 0.0;
        const auto ruleScvf = Dune::QuadratureRules<double, GG::GridView::dimension-1>::rule(scvfGeometry.type(), 3);
        for (const auto& qp : ruleScvf)
            volScvf += qp.weight()*GG::Extrusion::integrationElement(scvfGeometry, qp.position());

        if (!Dune::FloatCmp::eq(volScvf, GG::Extrusion::area(scvf)))
            DUNE_THROW(Dune::Exception, "Integration not correct! Integrated: " << volScvf << ", direct: " << GG::Extrusion::area(scvf));
    }

    // compare total volume/surface area to reference
    if (!Dune::FloatCmp::eq(refVolume, volume))
        DUNE_THROW(Dune::Exception, "Volume not correct! Reference: " << refVolume << ", computed: " << volume);
    if (!Dune::FloatCmp::eq(refSurface, surface))
        DUNE_THROW(Dune::Exception, "Surface not correct! Reference: " << refSurface << ", computed: " << surface);
}

} // end namespace Dumux

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // test disc extrusion 1d->2d
    {
        using Grid = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;

        struct GGTraits : public CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>
        { using Extrusion = RotationalExtrusion<0>; };
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

        std::cout << "Successfully tested disc extrusion." << std::endl;

    } // end disc extrusion

    // test spherical extrusion 1d->3d
    {
        using Grid = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
        struct GGTraits : public CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>
        { using Extrusion = SphericalExtrusion; };
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

        // compute the ball volume and the surface
        const double refVolume = 4.0/3.0*M_PI*(outerRadius*outerRadius*outerRadius - innerRadius*innerRadius*innerRadius);
        const double refSurface = 4.0*M_PI*(innerRadius*innerRadius + outerRadius*outerRadius);
        runTest(gg, refVolume, refSurface);

        std::cout << "Successfully tested spherical extrusion." << std::endl;

    } // end spherical extrusion

    // test rotational extrusion 2d->3d
    {
        // make a grid
        using Grid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
        using GlobalPosition = typename Grid::template Codim<0>::Geometry::GlobalCoordinate;
        const double innerRadius = 0.1;
        const double outerRadius = 1.0;
        const double height = 0.5;
        GlobalPosition lower({innerRadius, 0});
        GlobalPosition upper({outerRadius, height});
        std::array<unsigned int, Grid::dimension> els({20, 20});
        std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);
        auto leafGridView = grid->leafGridView();

        {
            // make GridGeometry
            struct GGTraits : public CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>
            { using Extrusion = RotationalExtrusion<0>; };
            using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, /*caching=*/false, GGTraits>;

            // make GridGeometry
            GridGeometry gg(leafGridView);
            gg.update();

            // compute the volume and the surface
            const auto centroidRadius = 0.5*(innerRadius + outerRadius);
            const auto width = (outerRadius - innerRadius);
            const double refVolume = height*width*2.0*M_PI*centroidRadius;
            const double refSurface = (2.0*height + 2.0*width)*2.0*M_PI*centroidRadius;
            runTest(gg, refVolume, refSurface);

            std::cout << "Successfully tested rotational extrusion with y-rotation axis." << std::endl;
        }

        {
            // make GridGeometry
            struct GGTraits : public CCTpfaDefaultGridGeometryTraits<typename Grid::LeafGridView>
            { using Extrusion = RotationalExtrusion<1>; };
            using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView, /*caching=*/false, GGTraits>;
            GridGeometry gg(leafGridView);
            gg.update();

            // compute the volume and the surface
            const auto centroidRadius = 0.5*height;
            const auto width = (outerRadius - innerRadius);
            const double refVolume = height*width*2.0*M_PI*centroidRadius;
            const double refSurface = 2.0*M_PI*height*width + 2*M_PI*height*height;
            runTest(gg, refVolume, refSurface);

            std::cout << "Successfully tested rotational extrusion with x-rotation axis." << std::endl;
        }

    } // end rotational extrusion

    return 0;
}
