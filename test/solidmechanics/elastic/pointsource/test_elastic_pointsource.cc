// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Test the point forces of the linear elastic momentum residual.
 *
 * Without displacement the residual of an element consists of its point forces only. A point force
 * \f$ \mathbf{q} \f$ at \f$ \mathbf{x} \f$ adds \f$ -\mathbf{q} \f$ to the control volume that contains
 * \f$ \mathbf{x} \f$ and \f$ -\mathbf{q} \varphi_i(\mathbf{x}) \f$ to a local dof \f$ i \f$ without control volume.
 * The basis functions \f$ \varphi_i \f$ of the PQ2 scheme are evaluated independently here.
 * A point on the face between two elements is split equally between them.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/pointsources.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/elastic/lameparams.hh>
#include <dumux/solidmechanics/elastic/model.hh>
#include <dumux/solidmechanics/elastic/spatialparams.hh>

namespace Dumux {

template<class Scalar, class GridGeometry>
class ElasticUnitSpatialParams
: public Experimental::ElasticSpatialParams<GridGeometry, Scalar, ElasticUnitSpatialParams<Scalar, GridGeometry>>
{
    using ParentType = Experimental::ElasticSpatialParams<GridGeometry, Scalar, ElasticUnitSpatialParams<Scalar, GridGeometry>>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using LameParams = Dumux::LameParams<Scalar>;

    ElasticUnitSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        lameParams_.setLambda(1.0);
        lameParams_.setMu(1.0);
    }

    const LameParams& lameParamsAtPos(const GlobalPosition&) const
    { return lameParams_; }

private:
    LameParams lameParams_;
};

template<class TypeTag>
class ElasticPointForceProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;

    struct Forces
    {
        std::vector<NumEqVector> values;

        template<class ElementDiscretization, class ElementVariables, class Context>
        NumEqVector operator()(const ElementDiscretization&, const ElementVariables&, const Context& context) const
        { return values[context.id()]; }
    };

    using Context = Experimental::IdPointSourceContext<Experimental::PointSourceContext<GlobalPosition>, std::size_t>;
    using PointSources = Experimental::PointSources<Context, Forces, LocalPosition>;

public:
    using ParentType::ParentType;

    void setPointSources(const std::vector<GlobalPosition>& positions, std::vector<NumEqVector> forces)
    {
        std::vector<Context> contexts;
        for (std::size_t i = 0; i < positions.size(); ++i)
            contexts.emplace_back(Experimental::PointSourceContext<GlobalPosition>(positions[i]), i);
        pointSources_ = Experimental::makePointSources(this->gridDiscretization(), contexts, Forces{std::move(forces)});
    }

    const PointSources& pointSources() const
    { return pointSources_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    { return BoundaryTypes(); }

    NumEqVector sourceAtPos(const GlobalPosition&) const
    { return NumEqVector(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

private:
    PointSources pointSources_;
};

namespace Properties {

namespace TTag {
struct ElasticPointForce { using InheritsFrom = std::tuple<Elastic>; };
struct ElasticPointForceBox { using InheritsFrom = std::tuple<ElasticPointForce, BoxModel>; };
struct ElasticPointForcePQ2Hybrid { using InheritsFrom = std::tuple<ElasticPointForce, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::ElasticPointForce> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::ElasticPointForce> { using type = ElasticPointForceProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ElasticPointForce>
{
    using type = ElasticUnitSpatialParams<GetPropType<TypeTag, Properties::Scalar>,
                                          GetPropType<TypeTag, Properties::GridGeometry>>;
};

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::ElasticPointForceBox> { static constexpr bool value = true; };

} // end namespace Properties
} // end namespace Dumux

namespace {

//! the Lagrange basis function of the given order with the given node, at xi, in the reference square
template<class LocalPosition>
double lagrange(int order, const LocalPosition& node, const LocalPosition& xi)
{
    double value = 1.0;
    for (int k = 0; k < 2; ++k)
    {
        const int j = std::lround(node[k]*order);
        for (int m = 0; m <= order; ++m)
            if (m != j)
                value *= (xi[k] - double(m)/order)/(double(j - m)/order);
    }
    return value;
}

template<class TypeTag>
double maxResidualDifference(int order)
{
    using namespace Dumux;
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
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(grid.leafGridView()))
    {
        elemDisc.bind(element);
        const auto elemVars = localView(gridVariables->curGridVars()).bind(element, elemDisc, x);
        const auto residual = localResidual.evalFluxAndSource(element, elemDisc, elemVars);

        std::vector<bool> isControlVolume(Detail::LocalDofs::numLocalDofs(elemDisc), false);
        for (const auto& scv : scvs(elemDisc))
            isControlVolume[scv.localDofIndex()] = true;

        const auto lower = element.geometry().corner(0);
        for (const auto& localDof : localDofs(elemDisc))
        {
            const auto node = ipData(elemDisc, localDof).local();
            GlobalPosition expected(0.0);
            for (std::size_t p = 0; p < positions.size(); ++p)
            {
                const auto relative = positions[p] - lower;
                if (!std::all_of(relative.begin(), relative.end(), [] (double r) { return r > -1e-12 && r < 0.25 + 1e-12; }))
                    continue;

                const bool onFace = std::any_of(relative.begin(), relative.end(), [] (double r) { return std::abs(r) < 1e-12 || std::abs(r - 0.25) < 1e-12; });
                auto xi = relative;
                xi /= 0.25;

                // the control volume of a corner is the part of the cell that is closer to it than to the other corners
                double weight = 1.0;
                if (isControlVolume[localDof.index()])
                    for (int k = 0; k < 2; ++k)
                        weight *= node[k] > 0.5 ? double(xi[k] > 0.5) : double(xi[k] < 0.5);
                else
                    weight = lagrange(order, node, xi);

                expected.axpy(-(onFace ? 0.5 : 1.0)*weight, forces[p]);
            }

            for (int k = 0; k < 2; ++k)
                maxDifference = std::max(maxDifference, std::abs(residual[localDof.index()][k] - expected[k]));
        }
    }
    return maxDifference;
}

} // end namespace

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    const auto box = maxResidualDifference<Properties::TTag::ElasticPointForceBox>(1);
    const auto hybrid = maxResidualDifference<Properties::TTag::ElasticPointForcePQ2Hybrid>(2);
    std::cout << "largest difference to the expected residual entries: box " << box << ", PQ2 hybrid " << hybrid << std::endl;

    if (box > 1e-13 || hybrid > 1e-13)
        DUNE_THROW(Dune::Exception, "The residual entries differ from the point force contributions");

    return 0;
}
