// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for rotation symmetric grid geometry
 */
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

template<class GG>
void runTest(const GG& gg, const double refVolume, const double refSurface)
{
    double volume = 0.0;
    double surface = 0.0;
    auto fvGeometry = localView(gg);
    for (const auto& element : elements(gg.gridView()))
    {
        fvGeometry.bind(element);

        for (const auto& scv : scvs(fvGeometry))
            volume += GG::Extrusion::volume(fvGeometry, scv);

        for (const auto& scvf : scvfs(fvGeometry))
            if (scvf.boundary())
                surface += GG::Extrusion::area(fvGeometry, scvf);

        // compare volume and integrated volume of one scv
        const auto& scv = *(scvs(fvGeometry).begin());
        const auto scvGeometry = fvGeometry.geometry(scv);

        double volScv = 0.0;
        const auto ruleScv = Dune::QuadratureRules<double, GG::GridView::dimension>::rule(scvGeometry.type(), 3);
        for (const auto& qp : ruleScv)
            volScv += qp.weight()*GG::Extrusion::integrationElement(scvGeometry, qp.position());

        if (!Dune::FloatCmp::eq(volScv, GG::Extrusion::volume(fvGeometry, scv)))
            DUNE_THROW(Dune::Exception, "Integration not correct! Integrated: " << volScv << ", direct: " << GG::Extrusion::volume(fvGeometry, scv));

        // compare area and integration area of one scvf
        const auto& scvf = *(scvfs(fvGeometry).begin());
        const auto scvfGeometry = fvGeometry.geometry(scvf);

        double volScvf = 0.0;
        const auto ruleScvf = Dune::QuadratureRules<double, GG::GridView::dimension-1>::rule(scvfGeometry.type(), 3);
        for (const auto& qp : ruleScvf)
            volScvf += qp.weight()*GG::Extrusion::integrationElement(scvfGeometry, qp.position());

        if (!Dune::FloatCmp::eq(volScvf, GG::Extrusion::area(fvGeometry, scvf)))
            DUNE_THROW(Dune::Exception, "Integration not correct! Integrated: " << volScvf << ", direct: " << GG::Extrusion::area(fvGeometry, scvf));
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

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

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
