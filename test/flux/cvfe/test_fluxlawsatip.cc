// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Test the flux laws evaluated at interpolation points with all control-volume finite element schemes
 *
 * All fields (temperature, densities, fractions, coefficients, displacement) are affine, so every
 * scheme represents them exactly, and the flux terms at every point of a sub-control volume face,
 * of an element and of a boundary face have to equal the exact flux terms:
 * - Fourier: \f$ -\lambda \nabla T \f$ with a scalar and a tensorial conductivity, interpolated
 *   or averaged over the face,
 * - Fick: \f$ -\varrho D^\kappa \nabla X^\kappa \f$ in both reference systems, without the main
 *   component and for a tracer fluid system,
 * - Maxwell-Stefan: the mass fluxes satisfy the Maxwell-Stefan relations of all components and sum to zero,
 * - Hooke: \f$ \lambda \operatorname{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2\mu \boldsymbol{\varepsilon} \f$
 *   for a displacement gradient that is not symmetric.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq1nonconforming.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>
#include <dumux/discretization/concepts.hh>
#include <dumux/discretization/defaultgridvariables.hh>
#include <dumux/discretization/cvfe/interpolate.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/flux/facetensoraverage.hh>
#include <dumux/flux/cvfe/coefficients_.hh>
#include <dumux/flux/cvfe/fickslaw_.hh>
#include <dumux/flux/cvfe/fourierslaw_.hh>
#include <dumux/flux/cvfe/hookeslaw_.hh>
#include <dumux/flux/cvfe/maxwellstefanslaw_.hh>
#include <dumux/solidmechanics/elastic/lameparams.hh>

namespace Dumux {

//! The fields carried by the solution
namespace Field {
enum : int {
    temperature,
    massDensity,
    molarDensity,
    massFraction,
    moleFraction = massFraction + 3,
    conductivity = moleFraction + 3,
    diffusionCoefficient,
    displacement,
    size = displacement + 3
};
} // end namespace Field

using FieldValues = Dune::FieldVector<double, Field::size>;
using Position = Dune::FieldVector<double, 3>;

//! The constant part and the gradient of every field
struct AffineFields
{
    AffineFields()
    {
        set(Field::temperature, 300.0, {2.0, -3.0, 1.5});
        set(Field::massDensity, 1000.0, {10.0, 5.0, -4.0});
        set(Field::molarDensity, 50.0, {1.0, -2.0, 0.5});
        set(Field::massFraction + 0, 0.2, {0.05, -0.02, 0.01});
        set(Field::massFraction + 1, 0.3, {-0.03, 0.04, -0.02});
        set(Field::moleFraction + 0, 0.25, {0.04, 0.03, -0.01});
        set(Field::moleFraction + 1, 0.35, {-0.05, 0.01, 0.02});
        set(Field::conductivity, 2.0, {0.5, 0.25, 0.1});
        set(Field::diffusionCoefficient, 1.5, {0.3, 0.2, 0.1});
        set(Field::displacement + 0, 0.01, {0.1, 0.2, -0.3});
        set(Field::displacement + 1, -0.02, {0.4, -0.1, 0.05});
        set(Field::displacement + 2, 0.03, {-0.2, 0.3, 0.1});

        // the fractions of the last component close the sum to one
        for (const int first : {Field::massFraction, Field::moleFraction})
        {
            constant[first + 2] = 1.0 - constant[first] - constant[first + 1];
            gradient[first + 2] = gradient[first];
            gradient[first + 2] += gradient[first + 1];
            gradient[first + 2] *= -1.0;
        }
    }

    template<class GlobalPosition>
    FieldValues operator()(const GlobalPosition& x) const
    {
        FieldValues values(constant);
        for (int i = 0; i < Field::size; ++i)
            for (int k = 0; k < GlobalPosition::dimension; ++k)
                values[i] += gradient[i][k]*x[k];
        return values;
    }

    template<class GlobalPosition>
    double value(int fieldIdx, const GlobalPosition& x) const
    { return (*this)(x)[fieldIdx]; }

    template<int dimWorld>
    Dune::FieldVector<double, dimWorld> grad(int fieldIdx) const
    {
        Dune::FieldVector<double, dimWorld> g;
        for (int k = 0; k < dimWorld; ++k)
            g[k] = gradient[fieldIdx][k];
        return g;
    }

    FieldValues constant;
    std::array<Position, Field::size> gradient;

private:
    void set(int fieldIdx, double c, const Position& g)
    { constant[fieldIdx] = c; gradient[fieldIdx] = g; }
};

const AffineFields fields;

struct TestFluidSystem
{
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 3;
    static constexpr bool isTracerFluidSystem() { return false; }
    static constexpr int getMainComponent(int phaseIdx) { return 0; }
    static double molarMass(int compIdx) { return std::array{0.018, 0.028, 0.044}[compIdx]; }
};

struct TestTracerFluidSystem
{
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 3;
    static constexpr bool isTracerFluidSystem() { return true; }
    static double molarMass(int compIdx) { return TestFluidSystem::molarMass(compIdx); }
};

//! The variables per local dof read the fields from the solution
class TestVariables
{
public:
    using PrimaryVariables = FieldValues;

    template<class ElementSolution, class Problem, class ElementDiscretization, class IpData>
    void update(const ElementSolution& elemSol, const Problem&, const ElementDiscretization&, const IpData& ipData)
    { values_ = elemSol[ipData.localDofIndex()]; }

    double temperature() const { return values_[Field::temperature]; }
    double density(int phaseIdx) const { return values_[Field::massDensity]; }
    double molarDensity(int phaseIdx) const { return values_[Field::molarDensity]; }
    double massFraction(int phaseIdx, int compIdx) const { return values_[Field::massFraction + compIdx]; }
    double moleFraction(int phaseIdx, int compIdx) const { return values_[Field::moleFraction + compIdx]; }
    double conductivity() const { return values_[Field::conductivity]; }
    double displacement(int dir) const { return values_[Field::displacement + dir]; }

    double averageMolarMass(int phaseIdx) const
    {
        double m = 0.0;
        for (int compIdx = 0; compIdx < TestFluidSystem::numComponents; ++compIdx)
            m += moleFraction(phaseIdx, compIdx)*TestFluidSystem::molarMass(compIdx);
        return m;
    }

    //! a coefficient that is affine in the fields and differs for every pair of components
    double effectiveDiffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    { return values_[Field::diffusionCoefficient]*(1.0 + compIIdx + compJIdx); }

private:
    FieldValues values_;
};

//! The exact diffusion coefficient of the test variables at a position
template<class GlobalPosition>
double diffusionCoefficient(const GlobalPosition& x, int compIIdx, int compJIdx)
{ return fields.value(Field::diffusionCoefficient, x)*(1.0 + compIIdx + compJIdx); }

struct TestSpatialParams
{
    template<class ElementDiscretization, class IpData>
    LameParams<double> lameParams(const ElementDiscretization&, const IpData&) const
    { return {3.0, 2.0}; }
};

struct TestProblem
{
    const TestSpatialParams& spatialParams() const { return spatialParams_; }
    std::string paramGroup() const { return ""; }
private:
    TestSpatialParams spatialParams_;
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
template<class G> struct FluxLawsTest {};
template<class G> struct TestBox { using InheritsFrom = std::tuple<FluxLawsTest<G>, BoxModel>; };
template<class G> struct TestPQ1FE { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ1FEModel>; };
template<class G> struct TestPQ1Bubble { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ1BubbleModel>; };
template<class G> struct TestPQ1BubbleHybrid { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ1BubbleHybridModel>; };
template<class G> struct TestPQ1BubbleFE { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ1BubbleFEModel>; };
template<class G> struct TestDiamond { using InheritsFrom = std::tuple<FluxLawsTest<G>, FaceCenteredDiamondModel>; };
template<class G> struct TestPQ1NonconformingFE { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ1NonconformingFEModel>; };
template<class G> struct TestPQ2Hybrid { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ2HybridModel>; };
template<class G> struct TestPQ2FE { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ2FEModel>; };
template<class G> struct TestPQ3Hybrid { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ3HybridModel>; };
template<class G> struct TestPQ3FE { using InheritsFrom = std::tuple<FluxLawsTest<G>, PQ3FEModel>; };
} // end namespace TTag

template<class TypeTag, class G>
struct Grid<TypeTag, TTag::FluxLawsTest<G>> { using type = G; };

template<class TypeTag, class G>
struct Scalar<TypeTag, TTag::FluxLawsTest<G>> { using type = double; };

template<class TypeTag, class G>
struct EnableGridGeometryCache<TypeTag, TTag::FluxLawsTest<G>> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

namespace Dumux {

struct ErrorTracker
{
    template<class Value, class Expected>
    void check(const Value& value, const Expected& expected, double scale, const char* what)
    {
        const double error = (value - expected).infinity_norm()/scale;
        maxError = std::max(maxError, error);
        if (!(error < tolerance))
            DUNE_THROW(Dune::Exception, what << ": relative error " << error << ", value " << value << ", expected " << expected);
    }

    static constexpr double tolerance = 1e-9;
    double maxError = 0.0;
};

//! Call f for each interpolation point at which the element variables of the scheme provide data
template<class ElementDiscretization, class F>
void forEachInterpolationPoint(const ElementDiscretization& elemDisc, F&& f)
{
    constexpr bool isFE = Experimental::Concepts::FEElementDiscretization<ElementDiscretization>;
    constexpr bool isHybrid = Experimental::Concepts::HybridElementDiscretization<ElementDiscretization>;

    if constexpr (!isFE)
        for (const auto& scvf : scvfs(elemDisc))
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
                f(qpData.ipData());

    if constexpr (isFE || isHybrid)
    {
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, elemDisc.element()))
            f(qpData.ipData());

        for (const auto& boundaryFace : boundaryFaces(elemDisc))
            for (const auto& qpData : CVFE::quadratureRule(elemDisc, boundaryFace))
                f(qpData.ipData());
    }
}

template<class TypeTag>
void checkFluxLaws(const typename GetPropType<TypeTag, Properties::Grid>::LeafGridView& gridView, const std::string& name)
{
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using Cache = Experimental::DefaultGridVariablesCache<GridDiscretization, TestProblem, TestVariables, false>;

    static constexpr int dim = GridDiscretization::GridView::dimension;
    static constexpr int dimWorld = GridDiscretization::GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<double, dimWorld>;
    using Tensor = Dune::FieldMatrix<double, dimWorld, dimWorld>;

    using FouriersLaw = Experimental::CVFEFouriersLawAtIp<double, GridDiscretization>;
    using FicksLawMass = Experimental::CVFEFicksLawAtIp<double, GridDiscretization, TestFluidSystem>;
    using FicksLawMolar = Experimental::CVFEFicksLawAtIp<double, GridDiscretization, TestFluidSystem,
                                                         ReferenceSystemFormulation::molarAveraged>;
    using FicksLawTracer = Experimental::CVFEFicksLawAtIp<double, GridDiscretization, TestTracerFluidSystem>;
    using MaxwellStefansLaw = Experimental::CVFEMaxwellStefansLawAtIp<double, GridDiscretization, TestFluidSystem>;
    using HookesLaw = Experimental::CVFEHookesLawAtIp<double, GridDiscretization>;

    constexpr int numComponents = TestFluidSystem::numComponents;
    constexpr int mainCompIdx = TestFluidSystem::getMainComponent(0);
    constexpr bool isFE = Experimental::Concepts::FEElementDiscretization<ElementDiscretization>;

    const GridDiscretization gridDisc(gridView);
    Dune::BlockVector<FieldValues> sol(gridDisc.numDofs());
    CVFE::InterpolationPolicy::L2Projection projection;
    projection.maxIterations = 1000;
    projection.residualReduction = 1e-16;
    CVFE::interpolate(gridDisc, sol, [](const auto& x) { return fields(x); }, projection);

    TestProblem problem;
    Cache gridCache(problem);
    gridCache.init(gridDisc, sol);

    // an anisotropic symmetric tensor
    Tensor anisotropy(0.0);
    for (int i = 0; i < dimWorld; ++i)
    {
        anisotropy[i][i] = 1.0 + i;
        if (i + 1 < dimWorld)
            anisotropy[i][i+1] = anisotropy[i+1][i] = 0.25;
    }

    const auto scalarConductivity = [] (const auto& vars) { return vars.conductivity(); };
    const auto tensorConductivity = [&] (const auto& vars) { auto t = anisotropy; t *= vars.conductivity(); return t; };
    const auto constantTensor = [&] (const auto&) { return anisotropy; };
    const auto gradT = fields.template grad<dimWorld>(Field::temperature);
    const double heatFluxScale = fields.constant[Field::conductivity]*gradT.infinity_norm();

    ErrorTracker errors;
    std::size_t numPoints = 0;
    auto elemDisc = localView(gridDisc);
    auto elemVars = localView(gridCache);
    for (const auto& element : elements(gridView))
    {
        elemDisc.bind(element);
        elemVars.bind(element, elemDisc, sol);

        forEachInterpolationPoint(elemDisc, [&](const auto& ipData)
        {
            ++numPoints;
            const auto& ipCache = cache(elemVars, ipData);
            const auto x = ipData.global();

            // Fourier
            {
                GlobalPosition expected = gradT;
                expected *= -fields.value(Field::conductivity, x);
                errors.check(FouriersLaw::fluxTerm(problem, element, elemDisc, elemVars, ipData, ipCache, scalarConductivity),
                             expected, heatFluxScale, "Fourier, scalar conductivity");

                anisotropy.mv(gradT, expected);
                expected *= -fields.value(Field::conductivity, x);
                errors.check(FouriersLaw::fluxTerm(problem, element, elemDisc, elemVars, ipData, tensorConductivity),
                             expected, heatFluxScale, "Fourier, tensorial conductivity");
            }

            if constexpr (!isFE)
            {
                if constexpr (Concept::ScvfIpData<std::decay_t<decltype(ipData)>>)
                {
                    GlobalPosition expected;
                    anisotropy.mv(gradT, expected);
                    expected *= -1.0;
                    errors.check(FouriersLaw::fluxTerm(problem, element, elemDisc, elemVars, ipData, ipCache,
                                                       Experimental::CVFE::faceTensorAverageAtScvf(constantTensor)),
                                 expected, heatFluxScale, "Fourier, face average of a constant tensor");

                    // the box dofs carry the field values at the vertices
                    if constexpr (GridDiscretization::discMethod == DiscretizationMethods::box)
                    {
                        const auto& scvf = elemDisc.scvf(ipData.scvfIndex());
                        const auto lambda = harmonicMean(
                            fields.value(Field::conductivity, elemDisc.scv(scvf.insideScvIdx()).dofPosition()),
                            fields.value(Field::conductivity, elemDisc.scv(scvf.outsideScvIdx()).dofPosition())
                        );
                        expected = gradT;
                        expected *= -lambda;
                        errors.check(FouriersLaw::fluxTerm(problem, element, elemDisc, elemVars, ipData, ipCache,
                                                           Experimental::CVFE::faceTensorAverageAtScvf(scalarConductivity)),
                                     expected, heatFluxScale, "Fourier, face average of the box values");
                    }
                }
            }

            // Fick
            const auto checkFick = [&]<class FicksLaw>(const FicksLaw&, int densityField, int fractionField, bool tracer)
            {
                const auto coefficient = [] (const auto& vars, int compIdx)
                { return FicksLaw::effectiveDiffusionCoefficient(vars, 0, compIdx); };
                const auto fluxTerms = FicksLaw::fluxTerm(problem, element, elemDisc, elemVars, 0, ipData, ipCache, coefficient);

                const double density = fields.value(densityField, x);
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    const auto gradX = fields.template grad<dimWorld>(fractionField + compIdx);
                    const double scale = density*diffusionCoefficient(x, 0, 0)*(1.0 + gradX.infinity_norm());
                    GlobalPosition expected(0.0);
                    if (tracer || compIdx != mainCompIdx)
                    {
                        expected = gradX;
                        expected *= -density*diffusionCoefficient(x, tracer ? 0 : mainCompIdx, compIdx);
                    }
                    errors.check(fluxTerms[compIdx], expected, scale, "Fick");
                }
            };
            checkFick(FicksLawMass{}, Field::massDensity, Field::massFraction, false);
            checkFick(FicksLawMolar{}, Field::molarDensity, Field::moleFraction, false);
            checkFick(FicksLawTracer{}, Field::massDensity, Field::massFraction, true);

            // Maxwell-Stefan: (rho/M) grad x_i = sum_j (x_i j_j/M_j - x_j j_i/M_i)/D_ij, sum_i j_i = 0
            {
                const auto coefficient = [] (const auto& vars, int i, int j)
                { return vars.effectiveDiffusionCoefficient(0, i, j); };
                const auto j = MaxwellStefansLaw::fluxTerm(problem, element, elemDisc, elemVars, 0, ipData, ipCache, coefficient);

                const double density = fields.value(Field::massDensity, x);
                double averageMolarMass = 0.0;
                for (int k = 0; k < numComponents; ++k)
                    averageMolarMass += fields.value(Field::moleFraction + k, x)*TestFluidSystem::molarMass(k);

                GlobalPosition sum(0.0);
                double fluxScale = 0.0;
                for (int i = 0; i < numComponents; ++i)
                {
                    sum += j[i];
                    fluxScale = std::max(fluxScale, j[i].infinity_norm());
                }
                errors.check(sum, GlobalPosition(0.0), fluxScale, "Maxwell-Stefan, sum of the fluxes");

                for (int i = 0; i < numComponents; ++i)
                {
                    const double xi = fields.value(Field::moleFraction + i, x);
                    GlobalPosition rhs(0.0);
                    for (int k = 0; k < numComponents; ++k)
                    {
                        if (k == i)
                            continue;
                        const double xk = fields.value(Field::moleFraction + k, x);
                        const double dik = diffusionCoefficient(x, i, k);
                        rhs.axpy(xi/(TestFluidSystem::molarMass(k)*dik), j[k]);
                        rhs.axpy(-xk/(TestFluidSystem::molarMass(i)*dik), j[i]);
                    }
                    auto lhs = fields.template grad<dimWorld>(Field::moleFraction + i);
                    lhs *= density/averageMolarMass;
                    errors.check(rhs, lhs, density/averageMolarMass*(1.0 + lhs.infinity_norm()), "Maxwell-Stefan relation");
                }
            }

            // Hooke
            if constexpr (dim == dimWorld)
            {
                Tensor gradU;
                for (int i = 0; i < dimWorld; ++i)
                    gradU[i] = fields.template grad<dimWorld>(Field::displacement + i);

                const auto lame = problem.spatialParams().lameParams(elemDisc, ipData);
                double trace = 0.0;
                for (int i = 0; i < dimWorld; ++i)
                    trace += gradU[i][i];

                Tensor expected(0.0);
                for (int i = 0; i < dimWorld; ++i)
                {
                    for (int k = 0; k < dimWorld; ++k)
                        expected[i][k] = lame.mu()*(gradU[i][k] + gradU[k][i]);
                    expected[i][i] += lame.lambda()*trace;
                }

                const auto sigma = HookesLaw::stressTensor(problem, element, elemDisc, elemVars, ipData, ipCache);
                for (int i = 0; i < dimWorld; ++i)
                    errors.check(sigma[i], expected[i], expected.infinity_norm(), "Hooke");
            }
        });
    }

    if (numPoints == 0)
        DUNE_THROW(Dune::Exception, name << ": no interpolation points checked");

    std::cout << name << ": " << numPoints << " interpolation points, maximum relative error " << errors.maxError << std::endl;
}

template<class Grid>
void checkAllSchemes(const Grid& grid, const std::string& gridName)
{
    using namespace Properties::TTag;
    const auto& gridView = grid.leafGridView();
    checkFluxLaws<TestBox<Grid>>(gridView, "box " + gridName);
    checkFluxLaws<TestPQ1Bubble<Grid>>(gridView, "pq1bubble " + gridName);
    checkFluxLaws<TestDiamond<Grid>>(gridView, "diamond " + gridName);
    checkFluxLaws<TestPQ1BubbleHybrid<Grid>>(gridView, "pq1bubble hybrid " + gridName);
    checkFluxLaws<TestPQ2Hybrid<Grid>>(gridView, "pq2 hybrid " + gridName);
    checkFluxLaws<TestPQ3Hybrid<Grid>>(gridView, "pq3 hybrid " + gridName);
    checkFluxLaws<TestPQ1FE<Grid>>(gridView, "pq1 fe " + gridName);
    checkFluxLaws<TestPQ1BubbleFE<Grid>>(gridView, "pq1bubble fe " + gridName);
    checkFluxLaws<TestPQ1NonconformingFE<Grid>>(gridView, "pq1nonconforming fe " + gridName);
    checkFluxLaws<TestPQ2FE<Grid>>(gridView, "pq2 fe " + gridName);
    checkFluxLaws<TestPQ3FE<Grid>>(gridView, "pq3 fe " + gridName);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    // anisotropic cells, so that the shape function gradients depend on the inverse Jacobian
    auto makeCubeGrid = []<int dim>(std::integral_constant<int, dim>)
    {
        using Grid = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim>>;
        std::array<std::vector<double>, dim> coordinates;
        for (int k = 0; k < dim; ++k)
            coordinates[k] = {0.0, 0.1*(k + 1), 0.5, 1.0 + k};
        return std::make_unique<Grid>(coordinates);
    };

    checkAllSchemes(*makeCubeGrid(std::integral_constant<int, 2>{}), "cube 2d");
    checkAllSchemes(*makeCubeGrid(std::integral_constant<int, 3>{}), "cube 3d");

#if HAVE_DUNE_ALUGRID
    auto makeSimplexGrid = []<int dim>(std::integral_constant<int, dim>)
    {
        using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
        Dune::FieldVector<double, dim> lower(0.0), upper;
        std::array<unsigned int, dim> cells;
        for (int k = 0; k < dim; ++k)
        {
            upper[k] = 1.0 + k;
            cells[k] = 2 + k;
        }
        return Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, cells);
    };

    checkAllSchemes(*makeSimplexGrid(std::integral_constant<int, 2>{}), "simplex 2d");
    checkAllSchemes(*makeSimplexGrid(std::integral_constant<int, 3>{}), "simplex 3d");
#endif

    return 0;
}
