// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief The inertia operator of the linear elastic model with the multi-stage assembler for
 *        local degrees of freedom, on the geometry of Cook's membrane
 *
 * For the box and the finite element schemes, the inertia operator of a translation summed over
 * all degrees of freedom is the density times the area times the translation. For all schemes,
 * the assembled residual is the weighted sum of the inertia, storage and spatial operators plus
 * the contribution of known states, and the assembled Jacobian agrees with central differences
 * of the residual. The problem is linear, so the numeric differentiation step can be chosen large
 * enough (Assembly.NumericDifference.BaseEpsilon) that its round-off does not limit the agreement.
 */
#include <config.h>

#include <iostream>
#include <memory>
#include <random>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/concepts.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/assembly/multistageassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>

#include <test/experimental/timestepping/implicitsolvecheck.hh>

#include "properties.hh"

namespace Dumux {

template<class GridGeometry, class Scalar>
class CooksMembraneInertiaSpatialParams
: public Experimental::ElasticSpatialParams<GridGeometry, Scalar, CooksMembraneInertiaSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = CooksMembraneInertiaSpatialParams<GridGeometry, Scalar>;
    using ParentType = Experimental::ElasticSpatialParams<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using LameParams = Dumux::LameParams<Scalar>;

    CooksMembraneInertiaSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , density_(getParam<Scalar>("SpatialParams.SolidDensity"))
    {
        const auto E = getParam<Scalar>("SpatialParams.YoungsModulus");
        const auto nu = getParam<Scalar>("SpatialParams.PoissonRatio");
        lameParams_.setLambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
        lameParams_.setMu(E/(2.0*(1.0 + nu)));
    }

    const LameParams& lameParamsAtPos(const GlobalPosition&) const
    { return lameParams_; }

    Scalar solidDensityAtPos(const GlobalPosition&) const
    { return density_; }

private:
    LameParams lameParams_;
    Scalar density_;
};

template<class TypeTag>
class CooksMembraneInertiaProblem : public CooksMembraneProblem<TypeTag>
{
public:
    using CooksMembraneProblem<TypeTag>::CooksMembraneProblem;

    static constexpr bool enableInertia()
    { return true; }
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct CooksMembraneInertia
{ using InheritsFrom = std::tuple<CooksMembrane>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembraneInertia>
{ using type = CooksMembraneInertiaProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembraneInertia>
{
    using type = CooksMembraneInertiaSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                                   GetPropType<TypeTag, Properties::Scalar>>;
};

} // end namespace Dumux::Properties

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::CooksMembraneInertia;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::MultiStageAssembler<TypeTag, DiffMethod::numeric>;
    static_assert(Assembler::hasInertia(), "The problem enables inertia");

    // the method only determines the Jacobian pattern here
    const auto method = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
    const auto xOld = x;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, method, xOld);

    bool passed = true;

    // a translation: for the box and the finite element schemes, whose basis functions sum to one
    // and whose control volumes or basis functions carry all of the mass, the inertia summed over
    // all degrees of freedom is rho |Omega| a
    using ElementDiscretization = typename GridGeometry::LocalView;
    constexpr bool isBox = std::is_same_v<typename GridGeometry::DiscretizationMethod, DiscretizationMethods::Box>;
    constexpr bool isFE = Experimental::Concepts::FEElementDiscretization<ElementDiscretization>;
    if constexpr (isBox || isFE)
    {
        const Dune::FieldVector<Scalar, 2> translation{0.3, -0.2};
        SolutionVector xTranslated(x.size());
        for (auto& block : xTranslated)
            block = translation;
        gridVariables->update(xTranslated);

        typename Assembler::ResidualType inertia;
        assembler->assembleOperators(xTranslated, Experimental::TimeLevel<Scalar>(0.0), {.inertia = &inertia});

        Dune::FieldVector<Scalar, 2> total(0.0);
        for (const auto& block : inertia)
            total += block;

        Scalar area = 0.0;
        for (const auto& element : elements(leafGridView))
            area += element.geometry().volume();

        auto expected = translation;
        expected *= problem->spatialParams().solidDensityAtPos({}) * area;
        auto difference = total;
        difference -= expected;
        std::cout << "[Inertia] translation: summed inertia " << total << ", expected " << expected << std::endl;
        if (!(difference.infinity_norm() < 1e-12*expected.infinity_norm()))
            passed = false;
    }

    // a random state
    std::mt19937 generator(42);
    std::uniform_real_distribution<Scalar> distribution(-1.0, 1.0);
    SolutionVector xRandom(x.size());
    for (auto& block : xRandom)
        for (auto& entry : block)
            entry = distribution(generator);
    gridVariables->update(xRandom);

    const auto numConstrained = Test::checkImplicitSolveContract(*assembler, xRandom, 0.0, 0.5, 1.0, 400.0);
    if (numConstrained == 0)
        passed = false;

    // the Jacobian of the inertia alone and of all three operators
    const Scalar inertiaDeviation = Test::jacobianDeviationFromFiniteDifferences(
        *assembler, xRandom, 0.0, Experimental::OperatorWeights<Scalar>{400.0, 0.0, 0.0}, 1e-6);
    const Scalar fullDeviation = Test::jacobianDeviationFromFiniteDifferences(
        *assembler, xRandom, 0.0, Experimental::OperatorWeights<Scalar>{400.0, 0.5, 1.0}, 1e-6);
    std::cout << "[Inertia] relative deviation of the Jacobian from central differences: inertia only "
              << inertiaDeviation << ", all operators " << fullDeviation << std::endl;
    if (!(inertiaDeviation < 1e-7) || !(fullDeviation < 1e-5))
        passed = false;

    std::cout << (passed ? "Inertia test: PASSED" : "Inertia test: FAILED") << std::endl;
    return passed ? 0 : 1;
}
