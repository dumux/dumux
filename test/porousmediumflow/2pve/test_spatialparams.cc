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

template<class GridGeometry, class SpatialParamsFine>
class SpatialParams
: public TwoPVESpatialParams<GridGeometry, double, SpatialParamsFine, SpatialParams<GridGeometry, SpatialParamsFine>>
{
    using ParentType = TwoPVESpatialParams<GridGeometry, double, SpatialParamsFine, SpatialParams<GridGeometry, SpatialParamsFine>>;
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
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SpatialParamsFine = TwoPVETest::FineSpatialParams<Element>;

    const Dune::FieldVector<Scalar, dim> lowerLeft({1.0, -2.0});
    const Dune::FieldVector<Scalar, dim> upperRight({5.0, 6.0});
    const std::array<int, dim> coarseCells({2, 1});
    const std::array<int, dim> fineCells({2, 4});
    const Scalar fineCellHeight = (upperRight[dim-1] - lowerLeft[dim-1])/fineCells[dim-1];

    Grid coarseGrid(lowerLeft, upperRight, coarseCells);
    Grid fineGrid(lowerLeft, upperRight, fineCells);
    auto coarseGridGeometry = std::make_shared<GridGeometry>(coarseGrid.leafGridView());
    auto fineGridGeometry = std::make_shared<GridGeometry>(fineGrid.leafGridView());
    const VEColumnMapping<GridGeometry, Scalar> columnMapping(coarseGridGeometry, fineGridGeometry);

    const TwoPVETest::SpatialParams<GridGeometry, SpatialParamsFine> spatialParams(
        coarseGridGeometry, columnMapping, std::make_shared<const SpatialParamsFine>(), fineCellHeight
    );

    const Scalar bottom = lowerLeft[dim-1];
    const Scalar top = upperRight[dim-1];
    const Scalar height = top - bottom;
    for (const auto& element : elements(coarseGridGeometry->gridView()))
    {
        // vertical averages of the fine-level fields over the column
        const Scalar x = element.geometry().center()[0];
        const Scalar expectedPermeability = 1.0e-12*((1.0 + x)*height + (top*top - bottom*bottom))/height;
        const Scalar expectedPorosity = ((0.1 + 0.01*x)*height + 0.01*(top*top - bottom*bottom))/height;

        TwoPVETest::checkClose(spatialParams.permeabilityAtElement(element), expectedPermeability, "column-averaged permeability");
        TwoPVETest::checkClose(spatialParams.porosityAtElement(element), expectedPorosity, "column-averaged porosity");

        const auto fvGeometry = localView(*coarseGridGeometry).bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const int elemSol = 0;
            TwoPVETest::checkClose(spatialParams.permeability(element, scv, elemSol), expectedPermeability, "scv permeability");
            TwoPVETest::checkClose(spatialParams.porosity(element, scv, elemSol), expectedPorosity, "scv porosity");
        }
    }

    if (spatialParams.template wettingPhaseAtPos<TwoPVETest::FluidSystem>(lowerLeft) != TwoPVETest::FluidSystem::phase0Idx)
        DUNE_THROW(Dune::Exception, "Expected the first phase to be the wetting phase");

    return 0;
}
