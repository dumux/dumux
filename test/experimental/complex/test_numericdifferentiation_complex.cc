// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Unit tests for numeric differentiation and the numeric epsilon with complex-valued variables
 */
#include <config.h>

#include <cmath>
#include <complex>
#include <iostream>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/ftraits.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/common/integrate.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>

namespace Dumux::Test {

template<class A, class B>
void checkClose(const A& a, const B& b, double tol, const std::string& what)
{
    using std::abs;
    const double diff = abs(a - b);
    if (diff > tol)
        DUNE_THROW(Dune::Exception, what << ": |" << a << " - " << b << "| = " << diff << " > " << tol);
}

} // end namespace Dumux::Test

int main(int argc, char** argv)
{
    using namespace Dumux;
    using Complex = std::complex<double>;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    // the epsilon has to be real-valued for complex-valued variables
    {
        const Complex x0(3.0, -4.0);
        const auto eps = NumericDifferentiation::epsilon(x0);
        static_assert(std::is_same_v<std::decay_t<decltype(eps)>, double>);
        Test::checkClose(eps, 1e-10*(5.0 + 1.0), 1e-20, "epsilon scaling with |x0|");

        const auto epsReal = NumericDifferentiation::epsilon(2.0);
        static_assert(std::is_same_v<std::decay_t<decltype(epsReal)>, double>);
        Test::checkClose(epsReal, 1e-10*3.0, 1e-20, "epsilon scaling for real x0");
    }

    // the numeric epsilon class parses a real base epsilon also for complex primary variables
    {
        const double baseEps = getParam<double>("Assembly.NumericDifference.BaseEpsilon", 1e-10);
        NumericEpsilon<Complex, 2> eps;
        const auto e = eps(Complex(1.0, 1.0), 0);
        static_assert(std::is_same_v<std::decay_t<decltype(e)>, double>);
        Test::checkClose(e, baseEps*(std::sqrt(2.0) + 1.0), 1e-12*baseEps, "NumericEpsilon for complex primary variable");

        NumericEpsilon<double, 2> epsReal;
        Test::checkClose(epsReal(2.0, 1), baseEps*3.0, 1e-12*baseEps, "NumericEpsilon for real primary variable");
    }

    // derivative of a holomorphic function with a real perturbation
    {
        const Complex a(2.0, 1.0), b(-1.0, 3.0);
        auto f = [&](Complex z){ return Dune::FieldVector<Complex, 1>(a*z*z + b*z); };
        const Complex z0(0.5, -0.25);
        const auto fz0 = f(z0);
        const Complex exact = 2.0*a*z0 + b;

        // one-sided differences are first order in eps, central and five-point ones cancel the
        // quadratic term exactly, so only roundoff remains for them
        const double eps = 1e-6;
        for (int method : {1, 0, -1, 5})
        {
            Dune::FieldVector<Complex, 1> derivative;
            NumericDifferentiation::partialDerivative(f, z0, derivative, fz0, eps, method);
            static_assert(std::is_same_v<std::decay_t<decltype(derivative[0])>, Complex>);
            Test::checkClose(derivative[0], exact, method == 1 || method == -1 ? 1e-5 : 1e-8,
                             "partial derivative with method " + std::to_string(method));
        }

        // default epsilon overload (forward differences, eps of 1e-10 is roundoff-limited)
        Dune::FieldVector<Complex, 1> derivative;
        NumericDifferentiation::partialDerivative(f, z0, derivative, fz0);
        Test::checkClose(derivative[0], exact, 1e-4, "partial derivative with default epsilon");
    }

    // the extrusion factor of the basic volume variables stays real
    {
        struct Traits { using PrimaryVariables = Dune::FieldVector<Complex, 2>; };
        using VolVars = BasicVolumeVariables<Traits>;
        static_assert(std::is_same_v<std::decay_t<decltype(std::declval<VolVars>().extrusionFactor())>, double>);
        static_assert(std::is_same_v<std::decay_t<decltype(std::declval<VolVars>().priVar(0))>, Complex>);

        struct RealTraits { using PrimaryVariables = Dune::FieldVector<double, 2>; };
        using RealVolVars = BasicVolumeVariables<RealTraits>;
        static_assert(std::is_same_v<std::decay_t<decltype(std::declval<RealVolVars>().extrusionFactor())>, double>);
    }

    // interpolation and gradients of a complex-valued box solution, and the real-valued L2 norm
    {
        using Grid = Dune::YaspGrid<2>;
        Grid grid({1.0, 1.0}, {4, 4});
        using GridGeometry = BoxFVGridGeometry<double, Grid::LeafGridView>;
        GridGeometry gridGeometry(grid.leafGridView());

        // linear field u = (1 + 2i) x + (3 - i) y, reproduced exactly by the box scheme
        using PrimaryVariables = Dune::FieldVector<Complex, 1>;
        using SolutionVector = Dune::BlockVector<PrimaryVariables>;
        const Complex ax(1.0, 2.0), ay(3.0, -1.0);
        SolutionVector sol(gridGeometry.numDofs());
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(grid.leafGridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                sol[scv.dofIndex()] = PrimaryVariables(ax*scv.dofPosition()[0] + ay*scv.dofPosition()[1]);
        }

        for (const auto& element : elements(grid.leafGridView()))
        {
            const auto geometry = element.geometry();
            const auto elemSol = elementSolution(element, sol, gridGeometry);
            const auto center = geometry.center();
            const auto value = evalSolution(element, geometry, gridGeometry, elemSol, center);
            static_assert(std::is_same_v<std::decay_t<decltype(value)>, PrimaryVariables>);
            Test::checkClose(value[0], ax*center[0] + ay*center[1], 1e-14, "evalSolution with complex values");

            const auto gradients = evalGradients(element, geometry, gridGeometry, elemSol, center);
            static_assert(std::is_same_v<std::decay_t<decltype(gradients[0][0])>, Complex>);
            Test::checkClose(gradients[0][0], ax, 1e-14, "evalGradients x component");
            Test::checkClose(gradients[0][1], ay, 1e-14, "evalGradients y component");
        }

        SolutionVector zero(sol.size());
        zero = 0.0;
        const auto norm = integrateL2Error(gridGeometry, sol, zero, 4);
        static_assert(std::is_same_v<std::decay_t<decltype(norm)>, double>);
        // |u|^2 = |ax|^2 x^2 + |ay|^2 y^2 + 2 Re(ax conj(ay)) x y integrated over the unit square
        const double expected = std::sqrt(std::norm(ax)/3.0 + std::norm(ay)/3.0 + 2.0*std::real(ax*std::conj(ay))/4.0);
        Test::checkClose(norm, expected, 1e-13, "complex L2 norm");
    }

    std::cout << "All checks passed" << std::endl;
    return 0;
}
