// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test the point forces of the finite element momentum residual.
 *
 * At rest, with a constant pressure, the residual of an element consists of its point forces only.
 * A point force \f$ \mathbf{q} \f$ at \f$ \mathbf{x} \f$ adds \f$ -\mathbf{q} \varphi_i(\mathbf{x}) \f$ to the
 * local dof \f$ i \f$ with the biquadratic basis function \f$ \varphi_i \f$, which is evaluated independently here.
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
#include <dumux/discretization/pq2.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>

namespace Dumux {

template<class TypeTag>
class RestProblem : public CVFENavierStokesMomentumProblem<TypeTag>
{
    using ParentType = CVFENavierStokesMomentumProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    struct Forces
    {
        std::vector<GlobalPosition> values;

        template<class ElementDiscretization, class ElementVariables, class Context>
        GlobalPosition operator()(const ElementDiscretization&, const ElementVariables&, const Context& context) const
        { return values[context.id()]; }
    };

    using Context = Experimental::IdPointSourceContext<Experimental::PointSourceContext<GlobalPosition>, std::size_t>;
    using PointSources = Experimental::PointSources<Context, Forces, LocalPosition>;

public:
    using ParentType::ParentType;

    void setPointSources(const std::vector<GlobalPosition>& positions, std::vector<GlobalPosition> forces)
    {
        std::vector<Context> contexts;
        for (std::size_t i = 0; i < positions.size(); ++i)
            contexts.emplace_back(Experimental::PointSourceContext<GlobalPosition>(positions[i]), i);
        pointSources_ = Experimental::makePointSources(this->gridDiscretization(), contexts, Forces{std::move(forces)});
    }

    const PointSources& pointSources() const
    { return pointSources_; }

    Scalar pressureAtPos(const GlobalPosition&) const
    { return 0.0; }

    Scalar densityAtPos(const GlobalPosition&) const
    { return 1.0; }

    Scalar effectiveViscosityAtPos(const GlobalPosition&) const
    { return 1.0; }

private:
    PointSources pointSources_;
};

namespace Properties {

namespace TTag {
struct RestMomentumFE
{
    using InheritsFrom = std::tuple<NavierStokesMomentumFE, PQ2FEModel>;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::RestMomentumFE>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::RestMomentumFE>
{ using type = RestProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RestMomentumFE>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

} // end namespace Properties
} // end namespace Dumux

//! the biquadratic Lagrange basis function with the given node, at xi, in the reference cube
template<class LocalPosition>
double biquadratic(const LocalPosition& node, const LocalPosition& xi)
{
    double value = 1.0;
    for (int k = 0; k < 2; ++k)
    {
        const int j = std::lround(node[k]*2);
        for (int m = 0; m <= 2; ++m)
            if (m != j)
                value *= (xi[k] - 0.5*m)/(0.5*(j - m));
    }
    return value;
}

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::RestMomentumFE;

    initialize(argc, argv);
    Parameters::init([](Dune::ParameterTree& params)
    {
        params["Problem.EnableGravity"] = "false";
        params["Problem.EnableInertiaTerms"] = "false";
        params["Component.LiquidDensity"] = "1";
        params["Component.LiquidKinematicViscosity"] = "1";
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
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using GlobalPosition = typename Grid::LeafGridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    const LocalResidual localResidual(problem.get());

    // a point in the interior of the cell [0.25, 0.5]^2 and a point on the face between this cell and its right neighbor
    const std::vector<GlobalPosition> positions{{0.31, 0.29}, {0.5, 0.4}};
    const std::vector<GlobalPosition> forces{{3.0, -1.0}, {5.0, 2.0}};
    problem->setPointSources(positions, forces);

    double maxDifference = 0.0;
    GlobalPosition total(0.0);
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(grid.leafGridView()))
    {
        elemDisc.bind(element);
        const auto elemVars = localView(gridVariables->curGridVars()).bind(element, elemDisc, x);
        const auto residual = localResidual.evalFluxAndSource(element, elemDisc, elemVars);

        const auto lower = element.geometry().corner(0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto node = ipData(elemDisc, localDof).local();
            GlobalPosition expected(0.0);
            for (std::size_t p = 0; p < positions.size(); ++p)
            {
                const auto relative = positions[p] - lower;
                if (std::all_of(relative.begin(), relative.end(), [] (double r) { return r > -1e-12 && r < 0.25 + 1e-12; }))
                {
                    const bool onFace = std::any_of(relative.begin(), relative.end(), [] (double r) { return std::abs(r) < 1e-12 || std::abs(r - 0.25) < 1e-12; });
                    auto xi = relative;
                    xi /= 0.25;
                    expected.axpy(-(onFace ? 0.5 : 1.0)*biquadratic(node, xi), forces[p]);
                }
            }

            for (int k = 0; k < 2; ++k)
            {
                maxDifference = std::max(maxDifference, std::abs(residual[localDof.index()][k] - expected[k]));
                total[k] += residual[localDof.index()][k];
            }
        }
    }

    std::cout << "largest difference to the expected residual entries: " << maxDifference << "\n"
              << "sum of all residual entries " << total << ", sum of the point forces " << forces[0] + forces[1] << std::endl;

    if (maxDifference > 1e-13)
        DUNE_THROW(Dune::Exception, "The residual entries differ from the point force contributions of the basis functions");

    if ((total + forces[0] + forces[1]).two_norm() > 1e-12)
        DUNE_THROW(Dune::Exception, "The point forces do not add up to their strengths");

    return 0;
}
