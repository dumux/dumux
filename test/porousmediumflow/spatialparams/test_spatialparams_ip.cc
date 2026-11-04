// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SpatialParameters
 * \brief Spatial parameters of porous-medium flow at interpolation points
 *
 * A parameter given as a function of the position has no default at an interpolation point, because a
 * local dof on an element face belongs to the elements on both sides. Calling one at an interpolation
 * point of spatial params that define it only by position fails to compile (`IPFUNCTION` 1 for
 * fluidMatrixInteraction, 2 wettingPhase, 3 characteristicLength, 4 factorEnergyTransfer,
 * 5 factorMassTransfer, 6 inertVolumeFraction). Constant parameters are unambiguous and defined there.
 */
#include <config.h>

#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/porousmediumflow/fvspatialparams1pconstant.hh>
#include <dumux/porousmediumflow/fvspatialparamsnonequilibrium.hh>

namespace Dumux::Test {

struct FluidSystem {};

struct TwoInertComponents
{
    static constexpr bool isInert() { return true; }
    static constexpr int numInertComponents = 2;
};

template<class GridGeometry, class Implementation>
class SpatialParamsByPosition
: public FVPorousMediumFlowSpatialParamsNonEquilibrium<GridGeometry, double, Implementation>
{
    using ParentType = FVPorousMediumFlowSpatialParamsNonEquilibrium<GridGeometry, double, Implementation>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = double;
    using ParentType::ParentType;

    double fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return 1.0 + globalPos[0]; }

    template<class FS>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return 1; }

    double characteristicLengthAtPos(const GlobalPosition& globalPos) const
    { return 2.0 + globalPos[0]; }

    double factorEnergyTransferAtPos(const GlobalPosition& globalPos) const
    { return 3.0 + globalPos[0]; }

    double factorMassTransferAtPos(const GlobalPosition& globalPos) const
    { return 4.0 + globalPos[0]; }

    template<class SolidSystem>
    double inertVolumeFractionAtPos(const GlobalPosition& globalPos, int compIdx) const
    { return 0.1*(compIdx + 1); }

    double porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.3; }

    double permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1e-12; }
};

//! Defines the parameters only by position
template<class GridGeometry>
class SpatialParamsAtPos
: public SpatialParamsByPosition<GridGeometry, SpatialParamsAtPos<GridGeometry>>
{
    using ParentType = SpatialParamsByPosition<GridGeometry, SpatialParamsAtPos<GridGeometry>>;
public:
    using ParentType::ParentType;
};

//! Defines the parameters at interpolation points in addition to the overloads of the base class
template<class GridGeometry>
class SpatialParamsAtIp
: public SpatialParamsByPosition<GridGeometry, SpatialParamsAtIp<GridGeometry>>
{
    using ParentType = SpatialParamsByPosition<GridGeometry, SpatialParamsAtIp<GridGeometry>>;
public:
    using ParentType::ParentType;
    using ParentType::fluidMatrixInteraction;
    using ParentType::wettingPhase;
    using ParentType::characteristicLength;
    using ParentType::factorEnergyTransfer;
    using ParentType::factorMassTransfer;

    template<class ElementDiscretization, class IpData, class ElementSolution>
    double fluidMatrixInteraction(const ElementDiscretization& elemDisc,
                                  const IpData& ipData,
                                  const ElementSolution& elemSol) const
    { return -this->fluidMatrixInteractionAtPos(ipData.global()); }

    template<class FS, class ElementDiscretization, class IpData, class ElementSolution>
    int wettingPhase(const ElementDiscretization& elemDisc,
                     const IpData& ipData,
                     const ElementSolution& elemSol) const
    { return 0; }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    double characteristicLength(const ElementDiscretization& elemDisc,
                                const IpData& ipData,
                                const ElementSolution& elemSol) const
    { return -this->characteristicLengthAtPos(ipData.global()); }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    double factorEnergyTransfer(const ElementDiscretization& elemDisc,
                                const IpData& ipData,
                                const ElementSolution& elemSol) const
    { return -this->factorEnergyTransferAtPos(ipData.global()); }

    template<class ElementDiscretization, class IpData, class ElementSolution>
    double factorMassTransfer(const ElementDiscretization& elemDisc,
                              const IpData& ipData,
                              const ElementSolution& elemSol) const
    { return -this->factorMassTransferAtPos(ipData.global()); }
};

void check(double value, double expected, const std::string& what)
{
    using std::abs;
    if (abs(value - expected) > 1e-14*abs(expected))
        DUNE_THROW(Dune::Exception, what << ": " << value << ", expected " << expected);
}

template<class GridGeometry>
void testSpatialParams(const GridGeometry& gridGeometry)
{
    const auto gg = std::make_shared<const GridGeometry>(gridGeometry);
    const SpatialParamsAtPos<GridGeometry> atPos(gg);
    const SpatialParamsAtIp<GridGeometry> atIp(gg);
    const FVPorousMediumFlowSpatialParamsOnePConstant<GridGeometry, double> constant(gg);
    const int elemSol = 0;

    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto center = scv.center();
            const auto dofPosition = scv.dofPosition();
            check(atPos.fluidMatrixInteraction(element, scv, elemSol), 1.0 + center[0], "fluidMatrixInteraction(scv)");
            check(atPos.template wettingPhase<FluidSystem>(element, scv, elemSol), 1, "wettingPhase(scv)");
            check(atPos.characteristicLength(element, scv, elemSol), 2.0 + dofPosition[0], "characteristicLength(scv)");
            check(atPos.factorEnergyTransfer(element, scv, elemSol), 3.0 + dofPosition[0], "factorEnergyTransfer(scv)");
            check(atPos.factorMassTransfer(element, scv, elemSol), 4.0 + dofPosition[0], "factorMassTransfer(scv)");
            check(atPos.template inertVolumeFraction<TwoInertComponents>(element, scv, elemSol, 1), 0.2, "inertVolumeFraction(scv)");
            check(atPos.porosity(element, scv, elemSol), 0.3, "porosity(scv)");
            check(atPos.permeability(element, scv, elemSol), 1e-12, "permeability(scv)");

            check(atIp.fluidMatrixInteraction(element, scv, elemSol), 1.0 + center[0], "fluidMatrixInteraction(scv), with ip overload");
            check(atIp.template wettingPhase<FluidSystem>(element, scv, elemSol), 1, "wettingPhase(scv), with ip overload");
            check(atIp.characteristicLength(element, scv, elemSol), 2.0 + dofPosition[0], "characteristicLength(scv), with ip overload");
            check(atIp.factorEnergyTransfer(element, scv, elemSol), 3.0 + dofPosition[0], "factorEnergyTransfer(scv), with ip overload");
            check(atIp.factorMassTransfer(element, scv, elemSol), 4.0 + dofPosition[0], "factorMassTransfer(scv), with ip overload");

            check(constant.porosity(element, scv, elemSol), 0.3, "constant porosity(scv)");
            check(constant.permeability(element, scv, elemSol), 2e-12, "constant permeability(scv)");
        }

        for (const auto& localDof : localDofs(fvGeometry))
        {
            const auto ip = ipData(fvGeometry, localDof);
            const auto globalPos = ip.global();

            check(atIp.fluidMatrixInteraction(fvGeometry, ip, elemSol), -1.0 - globalPos[0], "fluidMatrixInteraction(ip)");
            check(atIp.template wettingPhase<FluidSystem>(fvGeometry, ip, elemSol), 0, "wettingPhase(ip)");
            check(atIp.characteristicLength(fvGeometry, ip, elemSol), -2.0 - globalPos[0], "characteristicLength(ip)");
            check(atIp.factorEnergyTransfer(fvGeometry, ip, elemSol), -3.0 - globalPos[0], "factorEnergyTransfer(ip)");
            check(atIp.factorMassTransfer(fvGeometry, ip, elemSol), -4.0 - globalPos[0], "factorMassTransfer(ip)");

            check(constant.porosity(fvGeometry, ip, elemSol), 0.3, "constant porosity(ip)");
            check(constant.permeability(fvGeometry, ip, elemSol), 2e-12, "constant permeability(ip)");
            check(constant.temperature(fvGeometry, ip, elemSol), 300.0, "constant temperature(ip)");

#if IPFUNCTION == 1
            atPos.fluidMatrixInteraction(fvGeometry, ip, elemSol);
#elif IPFUNCTION == 2
            atPos.template wettingPhase<FluidSystem>(fvGeometry, ip, elemSol);
#elif IPFUNCTION == 3
            atPos.characteristicLength(fvGeometry, ip, elemSol);
#elif IPFUNCTION == 4
            atPos.factorEnergyTransfer(fvGeometry, ip, elemSol);
#elif IPFUNCTION == 5
            atPos.factorMassTransfer(fvGeometry, ip, elemSol);
#elif IPFUNCTION == 6
            atPos.template inertVolumeFraction<TwoInertComponents>(fvGeometry, ip, elemSol, 1);
#endif
        }
    }
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;
    initialize(argc, argv);
    Parameters::init([](Dune::ParameterTree& params)
    {
        params["Problem.EnableGravity"] = "false";
        params["SpatialParams.Porosity"] = "0.3";
        params["SpatialParams.Permeability"] = "2e-12";
        params["SpatialParams.Temperature"] = "300";
    });

    using Grid = Dune::YaspGrid<2>;
    using GridView = typename Grid::LeafGridView;
    Grid grid({1.0, 1.0}, {3, 2});
    Test::testSpatialParams(BoxFVGridGeometry<double, GridView, true>(grid.leafGridView()));

    std::cout << "Spatial parameters at interpolation points as expected" << std::endl;
    return 0;
}
