// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Placing gauging stations on a channel network and reading their discharge.
 *
 * A gauge is a surveyed point and the network a delineated one, so the two never coincide and
 * the station has to be snapped onto a face. Two ways of getting that wrong are silent: a
 * station can land on the wrong branch, which turns a tributary hydrograph into a copy of its
 * neighbour's, and it can be read with the sign of whichever side the mesh happens to call
 * inside, which makes the record incomparable with a gauge. Both are checked here, on a network
 * whose geometry is known rather than delineated, together with the discharge each station
 * reads in uniform flow.
 */
#include <config.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/discretization/box.hh>

#include <dumux/freeflow/shallowwater/longwave/model.hh>
#include <dumux/freeflow/shallowwater/longwave/gauges.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief A channel network with a bed falling along x and a uniform width and roughness
 */
template<class TypeTag>
class ChannelNetworkProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    struct SpatialParams
    {
        Scalar bedSurface(const Element&, const SubControlVolume& scv) const
        { return 100.0 - 0.1*scv.dofPosition()[0]; }

        Scalar manningN(const Element&) const
        { return 0.05; }

        template<class ElementSolution>
        Scalar extrusionFactor(const Element&, const SubControlVolume&, const ElementSolution&) const
        { return 2.0; }
    };

public:
    using ParentType::ParentType;

    const SpatialParams& spatialParams() const
    { return spatialParams_; }

private:
    SpatialParams spatialParams_;
};

namespace Properties::TTag {
struct ChannelNetwork
{
    using InheritsFrom = std::tuple<LongWave, BoxModel>;
    using Grid = Dune::FoamGrid<1, 2>;

    template<class TypeTag>
    using Problem = ChannelNetworkProblem<TypeTag>;
};
} // end namespace Properties::TTag

} // end namespace Dumux

namespace {

std::vector<std::string> failures;

void check(const std::string& name, bool condition, const std::string& detail = {})
{
    if (!condition)
        failures.push_back(name + (detail.empty() ? "" : ": " + detail));
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    Dune::MPIHelper::instance(argc, argv);
    Parameters::init([](Dune::ParameterTree& params){
        params["LongWave.Approximation"] = "diffusive";
    });

    using TypeTag = Properties::TTag::ChannelNetwork;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalPosition = Dune::FieldVector<double, 2>;

    // A Y: two branches meeting at (100, 0), then a trunk running east. The bed falls with x,
    // so downstream is +x on the trunk and towards the junction on both branches. Vertices are
    // inserted so that the two branches are ordered oppositely along the flow, which is what
    // makes the sign convention worth testing: taking the flux as positive out of the inside
    // sub-control volume would report one branch backwards.
    Dune::GridFactory<Grid> factory;
    factory.insertVertex({0.0, 50.0}); // 0, head of the north branch
    factory.insertVertex({50.0, 25.0}); // 1
    factory.insertVertex({100.0, 0.0}); // 2, the junction
    factory.insertVertex({50.0, -25.0}); // 3
    factory.insertVertex({0.0, -50.0}); // 4, head of the south branch
    factory.insertVertex({150.0, 0.0}); // 5
    factory.insertVertex({200.0, 0.0}); // 6, the outlet
    factory.insertElement(Dune::GeometryTypes::line, {0, 1});
    factory.insertElement(Dune::GeometryTypes::line, {1, 2});
    factory.insertElement(Dune::GeometryTypes::line, {3, 4}); // reversed on purpose
    factory.insertElement(Dune::GeometryTypes::line, {2, 3}); // reversed on purpose
    factory.insertElement(Dune::GeometryTypes::line, {2, 5});
    factory.insertElement(Dune::GeometryTypes::line, {5, 6});
    auto grid = factory.createGrid();

    auto gridGeometry = std::make_shared<GridGeometry>(grid->leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    // bed elevation per dof, as the spatial parameters define it
    std::vector<double> z(gridGeometry->numDofs(), 0.0);
    {
        auto fvGeometry = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                z[scv.dofIndex()] = problem->spatialParams().bedSurface(element, scv);
        }
    }

    using Gauges = LongWave::StreamGauges<GridGeometry>;
    const std::vector<std::string> names{"north", "south", "trunk"};
    const std::vector<GlobalPosition> positions{{25.0, 39.0}, {25.0, -39.0}, {175.0, 3.0}};

    // Each element of this network carries one interior face, at its midpoint. A station
    // surveyed a little off the channel must land on the midpoint of the reach it is beside.
    {
        const Gauges gauges(*gridGeometry, names, positions, z, 10.0);

        check("three stations located", gauges.size() == 3);

        std::vector<std::size_t> located;
        for (const auto& station : gauges.stations())
            located.push_back(station.element);
        check("the two branch stations are on different elements",
              located[0] != located[1],
              "both landed on element " + std::to_string(located[0]));
        check("the trunk station is on neither branch",
              located[2] != located[0] && located[2] != located[1]);

        // The north branch head is at (0, 50) and its first reach runs to (50, 25), so its
        // midpoint is (25, 37.5) -- 1.5 m from the surveyed (25, 39).
        check("north snapped to its own reach midpoint",
              std::abs(gauges.stations()[0].offset - 1.5) < 1e-9,
              "offset " + std::to_string(gauges.stations()[0].offset));

        // Both branches fall towards the junction, so both carry water the same way, but they
        // were inserted with opposite vertex order: the north element runs downhill from its
        // first vertex to its second, the south one uphill. The sign that makes both report
        // positive downstream is therefore +1 on one and -1 on the other.
        check("north falls from inside to outside, so downstream is out of the face",
              gauges.stations()[0].sign == 1,
              "sign " + std::to_string(gauges.stations()[0].sign));
        check("the reversed branch reports downstream with the opposite sign",
              gauges.stations()[1].sign == -1,
              "sign " + std::to_string(gauges.stations()[1].sign));
        check("the trunk falls away from the junction",
              gauges.stations()[2].sign == 1,
              "sign " + std::to_string(gauges.stations()[2].sign));

        const auto columns = gauges.columns();
        check("columns are named after the stations",
              columns.size() == 3 && columns[0] == "north[m^3/s]"
              && columns[2] == "trunk[m^3/s]");

        // At uniform depth the free surface parallels the bed, so every reach carries the
        // normal-flow discharge of Manning's law, b h^(5/3) sqrt(S) / n, downstream. The bed
        // falls by 0.1 per unit length in x, so the slope along a reach is 0.1 |dx| / length.
        const double depth = 0.2, width = 2.0, manningN = 0.05;
        SolutionVector sol(gridGeometry->numDofs());
        sol = depth;
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(sol);

        const auto discharges = gauges.discharges(*problem, *gridVariables, sol);
        const auto branchSlope = 0.1*50.0/std::hypot(50.0, 25.0);
        const std::vector<double> slopes{branchSlope, branchSlope, 0.1};
        for (std::size_t k = 0; k < names.size(); ++k)
        {
            const auto expected = width*std::pow(depth, 5.0/3.0)*std::sqrt(slopes[k])/manningN;
            check(names[k] + " reads the normal-flow discharge, positive downstream",
                  std::abs(discharges[k] - expected) < 1e-12*expected,
                  "got " + std::to_string(discharges[k]) + ", expected " + std::to_string(expected));
        }
    }

    // An empty network is not an error: a run that names no gauges simply has none.
    {
        const Gauges gauges(*gridGeometry, {}, {}, z, 10.0);
        check("no stations means empty", gauges.empty() && gauges.columns().empty());
    }

    // A station too far from any channel is an error rather than a snap. This is the case that
    // catches a survey given in the wrong coordinate system.
    {
        bool threw = false;
        try {
            const Gauges gauges(*gridGeometry, {"stray"}, {{25.0, 39.0}}, z, 1.0);
        } catch (const Dune::Exception&) { threw = true; }
        check("a station beyond the snap radius throws", threw);
    }

    // Mismatched inputs are caught rather than read past.
    {
        bool threw = false;
        try {
            const Gauges gauges(*gridGeometry, {"a", "b"}, {{25.0, 39.0}}, z, 10.0);
        } catch (const Dune::Exception&) { threw = true; }
        check("more names than positions throws", threw);

        threw = false;
        try {
            const Gauges gauges(*gridGeometry, {"a"}, {{25.0, 39.0}}, {0.0}, 10.0);
        } catch (const Dune::Exception&) { threw = true; }
        check("an elevation per dof is required", threw);
    }

    if (!failures.empty())
    {
        for (const auto& failure : failures)
            std::cerr << "FAILED " << failure << std::endl;
        return 1;
    }
    std::cout << "test_longwave_gauges: all checks passed" << std::endl;
    return 0;
}
