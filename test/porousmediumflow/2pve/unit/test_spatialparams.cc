// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Unit test for the coarse-level spatial parameters of the VE model.
 */

#include <config.h>

#include <array>
#include <cmath>
#include <memory>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>
#include <dumux/porousmediumflow/2pve/spatialparams.hh>

namespace Dumux::TwoPVETest {

// permeability and porosity are linear in both coordinates, so the midpoint rule integrates them exactly
template<class Element>
class FineSpatialParams
{
public:
    double permeabilityAtElement(const Element& element) const
    {
        const auto center = element.geometry().center();
        return 1.0e-12*(1.0 + center[0] + 2.0*center[1]);
    }

    double porosityAtElement(const Element& element) const
    {
        const auto center = element.geometry().center();
        return 0.1 + 0.01*center[0] + 0.02*center[1];
    }
};

template<class GridGeometry, class FineGridGeometry, class SpatialParamsFine>
class SpatialParams
: public TwoPVESpatialParams<GridGeometry, FineGridGeometry, double, SpatialParamsFine, SpatialParams<GridGeometry, FineGridGeometry, SpatialParamsFine>>
{
    using ParentType = TwoPVESpatialParams<GridGeometry, FineGridGeometry, double, SpatialParamsFine, SpatialParams<GridGeometry, FineGridGeometry, SpatialParamsFine>>;
public:
    using ParentType::ParentType;
};

struct FluidSystem
{
    static constexpr int phase0Idx = 0;
    static constexpr int phase1Idx = 1;
};

void checkClose(const double actual, const double expected, const std::string& quantity)
{
    using std::abs;
    if (abs(actual - expected) > 1.0e-14*abs(expected))
        DUNE_THROW(Dune::Exception, "Unexpected " << quantity << ": expected " << expected << ", obtained " << actual);
}

} // end namespace Dumux::TwoPVETest

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);
    Parameters::init([](auto& params){ params["Problem.EnableGravity"] = "true"; });

    static constexpr int dim = 2;
    using Scalar = double;
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<Scalar, dim>>;
    using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView>;
    using CoarseGrid = Dune::YaspGrid<dim-1, Dune::EquidistantOffsetCoordinates<Scalar, dim-1>>;
    using CoarseGridGeometry = CCTpfaFVGridGeometry<typename CoarseGrid::LeafGridView>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SpatialParamsFine = TwoPVETest::FineSpatialParams<Element>;

    const Dune::FieldVector<Scalar, dim> lowerLeft({1.0, -2.0});
    const Dune::FieldVector<Scalar, dim> upperRight({5.0, 6.0});
    const std::array<int, dim-1> coarseCells({2});
    const std::array<int, dim> fineCells({2, 4});
    const Scalar fineCellHeight = (upperRight[dim-1] - lowerLeft[dim-1])/fineCells[dim-1];

    CoarseGrid coarseGrid(Dune::FieldVector<Scalar, dim-1>({lowerLeft[0]}), Dune::FieldVector<Scalar, dim-1>({upperRight[0]}), coarseCells);
    Grid fineGrid(lowerLeft, upperRight, fineCells);
    auto coarseGridGeometry = std::make_shared<CoarseGridGeometry>(coarseGrid.leafGridView());
    auto fineGridGeometry = std::make_shared<GridGeometry>(fineGrid.leafGridView());
    const TwoPVEColumnMapping<CoarseGridGeometry, GridGeometry, Scalar> columnMapping(coarseGridGeometry, fineGridGeometry);

    const TwoPVETest::SpatialParams<CoarseGridGeometry, GridGeometry, SpatialParamsFine> spatialParams(
        coarseGridGeometry, columnMapping, std::make_shared<const SpatialParamsFine>(), fineCellHeight
    );

    const Scalar bottom = lowerLeft[dim-1];
    const Scalar top = upperRight[dim-1];
    const Scalar height = top - bottom;
    for (const auto& element : elements(coarseGridGeometry->gridView()))
    {
        // Vertical integrals over the fine-level fields.
        const Scalar x = element.geometry().center()[0];
        const Scalar expectedPermeabilityIntegral = 1.0e-12*((1.0 + x)*height + (top*top - bottom*bottom));
        const Scalar expectedPorosity = ((0.1 + 0.01*x)*height + 0.01*(top*top - bottom*bottom))/height;


        TwoPVETest::checkClose(spatialParams.permeabilityIntegralAtElement(element), expectedPermeabilityIntegral, "permeability integral");
        TwoPVETest::checkClose(spatialParams.porosityIntegralAtElement(element), expectedPorosity*height, "porosity integral");
        TwoPVETest::checkClose(spatialParams.columnHeightAtElement(element), height, "column height");
        if (spatialParams.gravity(element.geometry().center()).two_norm() != 0.0)
            DUNE_THROW(Dune::Exception, "Coarse horizontal gravity must vanish");
        if (std::abs(element.geometry().volume() - 2.0) > 1e-14)
            DUNE_THROW(Dune::Exception, "Coarse measure must be a length independent of height");

        const auto fvGeometry = localView(*coarseGridGeometry).bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const int elemSol = 0;
            TwoPVETest::checkClose(spatialParams.extrusionFactor(element, scv, elemSol), height, "extrusion height");
            TwoPVETest::checkClose(spatialParams.permeability(element, scv, elemSol), expectedPermeabilityIntegral/height, "scv average permeability");
            TwoPVETest::checkClose(spatialParams.porosity(element, scv, elemSol), expectedPorosity, "scv average porosity");
        }
    }

    if (spatialParams.template wettingPhaseAtPos<TwoPVETest::FluidSystem>(coarseGridGeometry->bBoxMin()) != TwoPVETest::FluidSystem::phase0Idx)
        DUNE_THROW(Dune::Exception, "Expected the first phase to be the wetting phase");

    return 0;
}
