// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Unit test for the grid and gravity requirements checked by the fine-level view of the VE model.
 */

#include <config.h>

#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/porousmediumflow/2p/indices.hh>
#include <dumux/porousmediumflow/2pve/finelevelview.hh>

namespace Dumux::TwoPVETest {

template<class GridGeometry>
class FineSpatialParams
{
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    FineSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry, GlobalPosition gravity = GlobalPosition({0.0, -9.81}))
    : gridGeometry_(gridGeometry), gravity_(gravity)
    {}

    double permeabilityAtElement(const Element& element) const
    { return 1.0e-12; }

    double porosityAtElement(const Element& element) const
    { return 0.2; }

    const GlobalPosition& gravity(const GlobalPosition&) const { return gravity_; }

    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    GlobalPosition gravity_;
};

template<class GlobalPosition>
class CoarseSpatialParams
{
    using MaterialLaw = FluidMatrix::BrooksCoreyNoReg<double>;
public:
    CoarseSpatialParams()
    : materialLaw_(typename MaterialLaw::BasicParams(1.0e5, 2.0), typename MaterialLaw::EffToAbsParams(0.1, 0.2))
    {}

    double temperatureAtPos(const GlobalPosition& globalPos) const
    { return 326.0; }

    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(materialLaw_); }

private:
    MaterialLaw materialLaw_;
};

template<class F>
void expectInvalidState(F&& function, const std::string& description)
{
    try
    {
        function();
    }
    catch (const Dune::InvalidStateException&)
    {
        return;
    }

    DUNE_THROW(Dune::Exception, "Expected an exception for " << description);
}

} // end namespace Dumux::TwoPVETest

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    static constexpr int dim = 2;
    using Grid = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim>>;
    using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView>;
    using GlobalPosition = Dune::FieldVector<double, dim-1>;
    using FinePosition = Dune::FieldVector<double, dim>;
    using CoarseGrid = Dune::YaspGrid<dim-1, Dune::TensorProductCoordinates<double, dim-1>>;
    using CoarseGridGeometry = CCTpfaFVGridGeometry<typename CoarseGrid::LeafGridView>;
    using FluidSystem = FluidSystems::TwoPImmiscible<double,
                                                     FluidSystems::OnePLiquid<double, Components::H2O<double>>,
                                                     FluidSystems::OnePGas<double, Components::CH4<double>>>;
    using SolutionVector = Dune::BlockVector<Dune::FieldVector<double, 2>>;
    using FineLevelView = TwoPVEFineLevelView<CoarseGridGeometry, GridGeometry, double, FluidSystem, TwoPIndices, SolutionVector, TwoPVETest::FineSpatialParams<GridGeometry>>;

    const std::vector<double> horizontalCoordinates({0.0, 1.0, 2.0});
    const auto makeFineLevelView = [&](const std::vector<double>& coarseCoordinates,
                                      const std::vector<double>& fineCoordinates,
                                      std::optional<FinePosition> gravity = std::nullopt)
    {
        auto coarseGrid = std::make_shared<CoarseGrid>(std::array<std::vector<double>, dim-1>({coarseCoordinates}));
        auto fineGrid = std::make_shared<Grid>(std::array<std::vector<double>, dim>({horizontalCoordinates, fineCoordinates}));
        auto coarseGridGeometry = std::make_shared<CoarseGridGeometry>(coarseGrid->leafGridView());
        auto fineGridGeometry = std::make_shared<GridGeometry>(fineGrid->leafGridView());
        using SpatialParams = TwoPVETest::FineSpatialParams<GridGeometry>;
        auto spatialParams = gravity
            ? std::make_shared<SpatialParams>(fineGridGeometry, *gravity)
            : std::make_shared<SpatialParams>(fineGridGeometry);
        return std::make_tuple(coarseGrid, fineGrid, coarseGridGeometry,
            std::make_shared<FineLevelView>(fineGridGeometry, coarseGridGeometry, spatialParams));
    };

    const auto [coarseGrid, fineGrid, coarseGridGeometry, fineLevelView] = makeFineLevelView(horizontalCoordinates, {0.0, 2.0, 4.0, 6.0});
    const auto coarseElement = *elements(coarseGridGeometry->gridView()).begin();
    const Dune::FieldVector<double, 2> waterSaturatedColumn({1.0e7, 0.0});
    const TwoPVETest::CoarseSpatialParams<GlobalPosition> coarseParams;
    const auto columnState = fineLevelView->makeColumnState(coarseElement, waterSaturatedColumn, coarseParams);
    if (std::abs(columnState.gravityNorm - 9.81) > 1.0e-14 || std::abs(columnState.gasPlumeDistance - 6.0) > 1.0e-14)
        DUNE_THROW(Dune::Exception, "Unexpected column state");

    const auto [inclinedCoarse, inclinedFine, inclinedGeometry, inclinedView] =
        makeFineLevelView(horizontalCoordinates, {0.0, 2.0, 4.0, 6.0}, FinePosition({1.0, -9.81}));
    TwoPVETest::expectInvalidState([&]{ inclinedView->makeColumnState(*elements(inclinedGeometry->gridView()).begin(), waterSaturatedColumn, coarseParams); },
                                   "gravity not aligned with the vertical axis");
    TwoPVETest::expectInvalidState([&]{ makeFineLevelView(horizontalCoordinates, {0.0, 2.0, 3.0, 6.0}); },
                                   "fine grid with non-uniform vertical spacing");
    TwoPVETest::expectInvalidState([&]{ makeFineLevelView({0.0, 2.0}, {0.0, 2.0, 4.0, 6.0}); },
                                   "coarse and fine grids with different horizontal subdivisions");
    return 0;
}
