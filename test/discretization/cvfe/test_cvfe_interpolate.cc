// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Tests for CVFE::interpolate: dof-position evaluation and L2 projection,
 *        for scalar- and vector-valued functions on cube and simplex grids.
 *
 * For each scheme and grid type, two functions are tested:
 *   - An in-space function (polynomial degree ≤ scheme order): both dof interpolation
 *     and L2 projection must reproduce it exactly (L2 error ≈ 0).
 *   - A function not in the approximation space (sin/cos/exp): the L2 projection
 *     error must be smaller than the dof interpolation error.
 *
 * Schemes and grids tested:
 *   - Box          – affine on a cube grid
 *   - Box          – affine on a simplex grid (requires dune-alugrid)
 *   - PQ1Bubble    – affine on a cube grid
 *   - PQ1Bubble    – affine on a simplex grid (requires dune-alugrid)
 *   - PQ2          – quadratic on a cube grid
 *   - PQ2          – quadratic on a simplex grid (requires dune-alugrid)
 *   - FCDiamond    – affine on a cube grid
 *   - FCDiamond    – affine on a simplex grid (requires dune-alugrid)
 *
 * Additionally, vector-valued (2-component) functions are tested on Box and PQ2
 * cube grids to verify that multi-component L2 projection works correctly.
 */
#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/io/format.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/facecentered/diamond/fvgridgeometry.hh>
#include <dumux/discretization/cvfe/interpolate.hh>

namespace Dumux {

/*!
 * \brief Compute the L2 error of a coefficient vector against an exact function.
 *
 * Works for both Scalar and vector-valued functions.
 */
template<class GridDiscretization, class Coeffs, class Function>
double computeL2Error(const GridDiscretization& gridDiscretization, const Coeffs& coeffs, Function&& fExact)
{
    static constexpr int dim = GridDiscretization::GridView::dimension;
    using Scalar = typename GridDiscretization::GridView::ctype;
    using ValueType = typename Coeffs::value_type;

    double l2error = 0.0;
    auto elemDisc = localView(gridDiscretization);
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        elemDisc.bind(element);
        const auto& elemGeo = elemDisc.elementGeometry();

        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
        const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(elemGeo.type(), 6);
        for (const auto& qp : quad)
        {
            elemDisc.feLocalBasis().evaluateFunction(qp.position(), shapeValues);

            ValueType interp(0.0);
            for (const auto& localDof : localDofs(elemDisc))
                interp += shapeValues[localDof.index()][0] * coeffs[localDof.dofIndex()];

            const auto diff = fExact(elemGeo.global(qp.position())) - interp;
            l2error += (diff * diff) * qp.weight() * elemGeo.integrationElement(qp.position());
        }
    }
    return std::sqrt(l2error);
}

/*!
 * \brief Run interpolation tests for a given grid geometry.
 *
 * Tests:
 * - Dof interpolation and L2 projection of fInSpace (must both be exact).
 * - Dof interpolation and L2 projection of fNotInSpace (L2 error of projection <= dof interpolation).
 */
template<class GridDiscretization, class FunctionInSpace, class FunctionNotInSpace>
void runTests(const GridDiscretization& gridDiscretization,
              FunctionInSpace&& fInSpace,
              FunctionNotInSpace&& fNotInSpace,
              const std::string& schemeName,
              const std::string& inSpaceName,
              const std::string& notInSpaceName,
              double toleranceExact = 1e-10)
{
    using Scalar = typename GridDiscretization::GridView::ctype;
    using GlobalPosition = Dune::FieldVector<Scalar, GridDiscretization::GridView::dimension>;
    using ValueTypeInSpace    = std::invoke_result_t<FunctionInSpace,    GlobalPosition>;
    using ValueTypeNotInSpace = std::invoke_result_t<FunctionNotInSpace, GlobalPosition>;

    // --- in-space function ---
    std::vector<ValueTypeInSpace> coeffsDofInSpace(gridDiscretization.numDofs(), ValueTypeInSpace(0.0));
    CVFE::interpolate(gridDiscretization, coeffsDofInSpace, fInSpace);

    std::vector<ValueTypeInSpace> coeffsL2InSpace(gridDiscretization.numDofs(), ValueTypeInSpace(0.0));
    CVFE::interpolate(gridDiscretization, coeffsL2InSpace, fInSpace, CVFE::InterpolationPolicy::L2Projection{});

    const double errDofInSpace = computeL2Error(gridDiscretization, coeffsDofInSpace, fInSpace);
    const double errL2InSpace  = computeL2Error(gridDiscretization, coeffsL2InSpace,  fInSpace);

    std::cout << Fmt::format("[{:<20}] dof-interp  {:<35}  L2 error = {:.3e}\n", schemeName, inSpaceName, errDofInSpace);
    std::cout << Fmt::format("[{:<20}] L2-proj     {:<35}  L2 error = {:.3e}\n", schemeName, inSpaceName, errL2InSpace);

    if (errDofInSpace > toleranceExact)
        DUNE_THROW(Dune::Exception, "[" << schemeName << "] Dof interpolation error for in-space function is not zero: " << errDofInSpace);
    if (errL2InSpace > toleranceExact)
        DUNE_THROW(Dune::Exception, "[" << schemeName << "] L2 projection error for in-space function is not zero: " << errL2InSpace);

    // --- not-in-space function ---
    std::vector<ValueTypeNotInSpace> coeffsDofNotInSpace(gridDiscretization.numDofs(), ValueTypeNotInSpace(0.0));
    CVFE::interpolate(gridDiscretization, coeffsDofNotInSpace, fNotInSpace);

    std::vector<ValueTypeNotInSpace> coeffsL2NotInSpace(gridDiscretization.numDofs(), ValueTypeNotInSpace(0.0));
    CVFE::interpolate(gridDiscretization, coeffsL2NotInSpace, fNotInSpace, CVFE::InterpolationPolicy::L2Projection{});

    const double errDofNotInSpace = computeL2Error(gridDiscretization, coeffsDofNotInSpace, fNotInSpace);
    const double errL2NotInSpace  = computeL2Error(gridDiscretization, coeffsL2NotInSpace,  fNotInSpace);

    std::cout << Fmt::format("[{:<20}] dof-interp  {:<35}  L2 error = {:.3e}\n", schemeName, notInSpaceName, errDofNotInSpace);
    std::cout << Fmt::format("[{:<20}] L2-proj     {:<35}  L2 error = {:.3e}\n", schemeName, notInSpaceName, errL2NotInSpace);

    if (errL2NotInSpace > errDofNotInSpace + toleranceExact)
        DUNE_THROW(Dune::Exception,
            "[" << schemeName << "] L2 projection error (" << errL2NotInSpace
            << ") is larger than dof interpolation error (" << errDofNotInSpace << ")");
}

} // namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    initialize(argc, argv);

    using Grid     = Dune::YaspGrid<3>;
    using GridView = Grid::LeafGridView;
    using Scalar   = double;
    using GlobalPosition = Dune::FieldVector<Scalar, 3>;

    // build a 2×2×2 unit-cube grid
    GlobalPosition lower(0.0), upper(1.0);
    std::array<unsigned int, 3> cells{{2, 2, 2}};
    auto grid         = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, cells);
    auto leafGridView = grid->leafGridView();

    // affine test function: f(x) = x[0] + 2·x[1] + 3·x[2]
    auto fAffine = [](const GlobalPosition& x) -> Scalar
    {
        Scalar val = 0.0;
        for (int i = 0; i < 3; ++i)
            val += Scalar(i + 1) * x[i];
        return val;
    };

    // quadratic test function: f(x) = x[0]² + 2·x[1]² + 3·x[2]²
    auto fQuadratic = [](const GlobalPosition& x) -> Scalar
    {
        Scalar val = 0.0;
        for (int i = 0; i < 3; ++i)
            val += Scalar(i + 1) * x[i] * x[i];
        return val;
    };

    // complex function not in any low-order approximation space
    auto fNotInSpace = [](const GlobalPosition& x) -> Scalar
    {
        using std::sin; using std::cos; using std::exp;
        constexpr Scalar pi = M_PI;
        return sin(pi*x[0]) * cos(pi*x[1]) * exp(x[2]);
    };

    {   using GG = BoxFVGridGeometry<Scalar, GridView, true>;
        GG gg(leafGridView);
        runTests(gg, fAffine, fNotInSpace, "Box-Cube", "x + 2y + 3z", "sin(px)cos(py)exp(z)"); }

    {   using GG = PQ1BubbleFVGridGeometry<Scalar, GridView, true>;
        GG gg(leafGridView);
        runTests(gg, fAffine, fNotInSpace, "PQ1Bubble-Cube", "x + 2y + 3z", "sin(px)cos(py)exp(z)"); }

    {   using GG = PQ2FVGridGeometry<Scalar, GridView, true>;
        GG gg(leafGridView);
        runTests(gg, fQuadratic, fNotInSpace, "PQ2-Cube", "x2 + 2y2 + 3z2", "sin(px)cos(py)exp(z)"); }

    {   using GG = FaceCenteredDiamondFVGridGeometry<GridView, true>;
        GG gg(leafGridView);
        runTests(gg, fAffine, fNotInSpace, "FCDiamond-Cube", "x + 2y + 3z", "sin(px)cos(py)exp(z)"); }

#if HAVE_DUNE_ALUGRID
    {
        using SimplexGrid     = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
        using SimplexGridView = SimplexGrid::LeafGridView;

        GlobalPosition lowerS(0.0), upperS(1.0);
        std::array<unsigned int, 3> cellsS{{2, 2, 2}};
        auto simplexGrid     = Dune::StructuredGridFactory<SimplexGrid>::createSimplexGrid(lowerS, upperS, cellsS);
        auto simplexGridView = simplexGrid->leafGridView();

        {   using GG = BoxFVGridGeometry<Scalar, SimplexGridView, true>;
            GG gg(simplexGridView);
            runTests(gg, fAffine, fNotInSpace, "Box-Simplex", "x + 2y + 3z", "sin(px)cos(py)exp(z)"); }

        {   using GG = PQ1BubbleFVGridGeometry<Scalar, SimplexGridView, true>;
            GG gg(simplexGridView);
            runTests(gg, fAffine, fNotInSpace, "PQ1Bubble-Simplex", "x + 2y + 3z", "sin(px)cos(py)exp(z)"); }

        {   using GG = PQ2FVGridGeometry<Scalar, SimplexGridView, true>;
            GG gg(simplexGridView);
            runTests(gg, fQuadratic, fNotInSpace, "PQ2-Simplex", "x2 + 2y2 + 3z2", "sin(px)cos(py)exp(z)"); }

        {   using GG = FaceCenteredDiamondFVGridGeometry<SimplexGridView, true>;
            GG gg(simplexGridView);
            runTests(gg, fAffine, fNotInSpace, "FCDiamond-Simplex", "x + 2y + 3z", "sin(px)cos(py)exp(z)"); }
    }
#endif

    // --- vector-valued tests (2 components) ---
    // in-space: both components affine
    auto fAffineVec = [](const GlobalPosition& x) -> Dune::FieldVector<Scalar, 2>
    { return {x[0] + 2*x[1] + 3*x[2], 2*x[0] + x[1] - x[2]}; };
    // not-in-space: trigonometric × exponential
    auto fNotInSpaceVec = [](const GlobalPosition& x) -> Dune::FieldVector<Scalar, 2>
    {
        using std::sin; using std::cos; using std::exp;
        constexpr Scalar pi = M_PI;
        return {sin(pi*x[0]) * cos(pi*x[1]) * exp(x[2]),
                cos(pi*x[0]) * sin(pi*x[1]) * exp(x[2])};
    };

    {   using GG = BoxFVGridGeometry<Scalar, GridView, true>;
        GG gg(leafGridView);
        runTests(gg, fAffineVec, fNotInSpaceVec,
            "Box-Cube (vec2)", "(x+2y+3z, 2x+y-z)", "(sin*cos*exp, cos*sin*exp)"); }

    {   using GG = PQ2FVGridGeometry<Scalar, GridView, true>;
        GG gg(leafGridView);
        auto fQuadraticVec = [](const GlobalPosition& x) -> Dune::FieldVector<Scalar, 2>
        {
            Scalar v0 = 0.0, v1 = 0.0;
            for (int i = 0; i < 3; ++i) { v0 += Scalar(i+1)*x[i]; v1 += Scalar(i+1)*x[i]*x[i]; }
            return {v0, v1};
        };
        runTests(gg, fQuadraticVec, fNotInSpaceVec,
            "PQ2-Cube (vec2)", "(x+2y+3z, x²+2y²+3z²)", "(sin*cos*exp, cos*sin*exp)"); }

    std::cout << "All tests passed.\n";
    return 0;
}
