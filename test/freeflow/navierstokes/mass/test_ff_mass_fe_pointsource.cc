// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test the point sources of the finite element continuity residual.
 *
 * A uniform velocity field and a constant density satisfy the continuity equation exactly,
 * so the residual of an element consists of its point sources only. A point source \f$ q \f$ at
 * \f$ \mathbf{x} \f$ adds \f$ -q \varphi_i(\mathbf{x}) \f$ to the local dof \f$ i \f$ with the
 * bilinear basis function \f$ \varphi_i \f$, which is evaluated independently here.
 * A point on the face between two elements is split equally between them.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/pointsources.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

template<class TypeTag>
class UniformFlowPointSourceProblem : public CVFENavierStokesMassProblem<TypeTag>
{
    using ParentType = CVFENavierStokesMassProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    struct Strengths
    {
        std::vector<Scalar> values;

        template<class ElementDiscretization, class ElementVariables, class Context>
        auto operator()(const ElementDiscretization&, const ElementVariables&, const Context& context) const
        { return Dune::FieldVector<Scalar, 1>(values[context.id()]); }
    };

    using Context = Experimental::IdPointSourceContext<Experimental::PointSourceContext<GlobalPosition>, std::size_t>;
    using PointSources = Experimental::PointSources<Context, Strengths, LocalPosition>;

public:
    using typename ParentType::BoundaryTypes;
    using ParentType::ParentType;

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    void setPointSources(const std::vector<GlobalPosition>& positions, std::vector<Scalar> strengths)
    {
        std::vector<Context> contexts;
        for (std::size_t i = 0; i < positions.size(); ++i)
            contexts.emplace_back(Experimental::PointSourceContext<GlobalPosition>(positions[i]), i);
        pointSources_ = Experimental::makePointSources(this->gridDiscretization(), contexts, Strengths{std::move(strengths)});
    }

    const PointSources& pointSources() const
    { return pointSources_; }

    GlobalPosition velocityAtPos(const GlobalPosition& globalPos) const
    { return {1.0, 0.5}; }

private:
    PointSources pointSources_;
};

namespace Properties {

namespace TTag {
struct UniformFlowPointSourceFE
{
    using InheritsFrom = std::tuple<NavierStokesMassOnePFE, PQ1FEModel>;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::UniformFlowPointSourceFE>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::UniformFlowPointSourceFE>
{ using type = UniformFlowPointSourceProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::UniformFlowPointSourceFE>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::UniformFlowPointSourceFE;

    initialize(argc, argv);
    Parameters::init([](Dune::ParameterTree& params)
    {
        params["Problem.EnableGravity"] = "false";
        params["Component.LiquidDensity"] = "1000";
        params["Component.LiquidKinematicViscosity"] = "1e-6";
        params["SpatialParams.Temperature"] = "293.15";
    });

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    const Grid grid({1.0, 1.0}, {4, 4});

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 1e5;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using GlobalPosition = typename Grid::LeafGridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    const LocalResidual localResidual(problem.get());

    // a point in the interior of the cell [0.25, 0.5]^2 and a point on the face between this cell and its right neighbor
    const std::vector<GlobalPosition> positions{{0.31, 0.29}, {0.5, 0.4}};
    const std::vector<double> strengths{3.0, 5.0};
    problem->setPointSources(positions, strengths);

    // the bilinear basis function of the vertex at the given corner of the cell with the given lower left corner
    const auto bilinear = [] (std::array<int, 2> corner, const GlobalPosition& lower, const GlobalPosition& position)
    {
        double value = 1.0;
        for (int k = 0; k < 2; ++k)
        {
            const double xi = (position[k] - lower[k])/0.25;
            value *= corner[k] ? xi : 1.0 - xi;
        }
        return value;
    };

    double maxDifference = 0.0, total = 0.0;
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(grid.leafGridView()))
    {
        elemDisc.bind(element);
        const auto elemVars = localView(gridVariables->curGridVars()).bind(element, elemDisc, x);
        const auto residual = localResidual.evalFluxAndSource(element, elemDisc, elemVars);

        const auto lower = element.geometry().corner(0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto corner = ipData(elemDisc, localDof).local();
            double expected = 0.0;
            for (std::size_t p = 0; p < positions.size(); ++p)
            {
                const auto relative = positions[p] - lower;
                // the point lies in the closure of the cell if all relative coordinates are in [0, 0.25]
                if (std::all_of(relative.begin(), relative.end(), [] (double r) { return r > -1e-12 && r < 0.25 + 1e-12; }))
                {
                    const bool onFace = std::any_of(relative.begin(), relative.end(), [] (double r) { return std::abs(r) < 1e-12 || std::abs(r - 0.25) < 1e-12; });
                    expected -= (onFace ? 0.5 : 1.0)*strengths[p]
                                *bilinear({int(corner[0] > 0.5), int(corner[1] > 0.5)}, lower, positions[p]);
                }
            }
            maxDifference = std::max(maxDifference, std::abs(residual[localDof.index()][0] - expected));
            total += residual[localDof.index()][0];
        }
    }

    std::cout << "largest difference to the expected residual entries: " << maxDifference << "\n"
              << "sum of all residual entries " << total << ", sum of the strengths " << std::accumulate(strengths.begin(), strengths.end(), 0.0) << std::endl;

    if (maxDifference > 1e-13)
        DUNE_THROW(Dune::Exception, "The residual entries differ from the point source contributions of the basis functions");

    if (std::abs(total + std::accumulate(strengths.begin(), strengths.end(), 0.0)) > 1e-12)
        DUNE_THROW(Dune::Exception, "The point sources do not add up to their strengths");

    return 0;
}
