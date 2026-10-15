// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief The inertia operator of the hyperelastic model with the multi-stage assembler for local
 *        degrees of freedom
 *
 * Checks that the inertia operator of a translation sums to the density times the volume times
 * the translation, that the assembled residual is the weighted sum of the inertia, storage and
 * spatial operators plus the contribution of known states, and that the assembled Jacobian agrees
 * with central differences of the residual. With CHECK_FIRST_ORDER_STEPPER_REJECTION, a
 * first-order time stepper is instantiated for the problem with inertia, which must not compile.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <random>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/solidmechanics/hyperelastic/spatialparams.hh>
#include <dumux/assembly/multistageassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>

#if CHECK_FIRST_ORDER_STEPPER_REJECTION
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#endif

#include <test/experimental/timestepping/implicitsolvecheck.hh>

#include "properties.hh"

namespace Dumux {

template<class TypeTag>
class HyperelasticityInertiaProblem : public HyperelasticityProblem<TypeTag>
{
    using ParentType = HyperelasticityProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using ParentType::ParentType;

    static constexpr bool enableInertia()
    { return true; }

    template<class ElementDiscretization, class IpData>
    Scalar solidDensity(const ElementDiscretization&, const IpData&) const
    { return this->spatialParams().solidDensity(); }
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct HyperelasticityInertiaTest
{
    using InheritsFrom = std::tuple<HyperelasticityTest>;

    template<class TypeTag>
    using Problem = HyperelasticityInertiaProblem<TypeTag>;
};
} // end namespace TTag

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HyperelasticityInertiaTest>
{
    using type = DefaultDynamicHyperelasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                                         GetPropType<TypeTag, Properties::Scalar>>;
};

} // end namespace Dumux::Properties

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::HyperelasticityInertiaTest;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

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

#if CHECK_FIRST_ORDER_STEPPER_REJECTION
    using LinearSolver = ILUBiCGSTABIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Newton = NewtonSolver<Assembler, LinearSolver>;
    auto newton = std::make_shared<Newton>(assembler, std::make_shared<LinearSolver>());
    Experimental::MultiStageTimeStepper<Newton> timeStepper(newton, method);
#endif

    bool passed = true;

    // a translation: the inertia summed over all degrees of freedom is rho |Omega| a
    {
        const Dune::FieldVector<Scalar, 3> translation{0.3, -0.2, 0.1};
        SolutionVector xTranslated(x.size());
        for (auto& block : xTranslated)
            block = translation;
        gridVariables->update(xTranslated);

        typename Assembler::ResidualType inertia;
        assembler->assembleOperators(xTranslated, Experimental::TimeLevel<Scalar>(0.0), {.inertia = &inertia});

        Dune::FieldVector<Scalar, 3> total(0.0);
        for (const auto& block : inertia)
            total += block;

        Scalar volume = 0.0;
        for (const auto& element : elements(leafGridView))
            volume += element.geometry().volume();

        auto expected = translation;
        expected *= problem->spatialParams().solidDensity()*volume;
        auto difference = total;
        difference -= expected;
        std::cout << "[Inertia] translation: summed inertia " << total << ", expected " << expected << std::endl;
        if (!(difference.infinity_norm() < 1e-12*expected.infinity_norm()))
            passed = false;
    }

    // a random state with nonzero stresses
    std::mt19937 generator(42);
    std::uniform_real_distribution<Scalar> distribution(-0.05, 0.05);
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
