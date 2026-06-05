//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/io/json.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>

namespace Dumux {

using Scalar = double;
using Complex = std::complex<Scalar>;
using Method = std::shared_ptr<const Experimental::MultiStageMethod<Scalar>>;
using ComplexFunction = std::function<Complex(Complex)>;

Complex stabilityFunction(const Method& method, const Complex z)
{
    std::vector<Complex> stage(method->numStages() + 1, 0.0);
    stage[0] = 1.0;

    for (std::size_t i = 1; i <= method->numStages(); ++i)
    {
        Complex rhs = 0.0;
        for (std::size_t k = 0; k < i; ++k)
            rhs += (method->temporalWeight(i, k) - z * method->spatialWeight(i, k)) * stage[k];

        const auto denominator = method->temporalWeight(i, i) - z * method->spatialWeight(i, i);
        stage[i] = -rhs / denominator;
    }

    return stage.back();
}

/*!
 * \brief Analytical stability function of an explicit Runge-Kutta method
 *        of classical order p with p stages: R(z) = sum_{k=0}^{p} z^k/k!
 */
Complex truncatedExponential(const Complex& z, const std::size_t order)
{
    Complex term = 1.0;
    Complex sum = 1.0;
    for (std::size_t k = 1; k <= order; ++k)
    {
        term *= z / static_cast<Scalar>(k);
        sum += term;
    }
    return sum;
}

void expect(const bool cond, const std::string& msg)
{
    if (!cond)
        DUNE_THROW(Dune::InvalidStateException, msg);
}

void expectNear(const Scalar value, const Scalar reference, const Scalar tolerance, const std::string& msg)
{
    using std::abs;
    if (abs(value - reference) > tolerance)
        DUNE_THROW(Dune::InvalidStateException, msg);
}

//! Central-difference derivative (R is holomorphic, so a real step direction suffices)
Complex derivative(const ComplexFunction& R, const Complex& z, const Scalar h = 1.0e-6)
{
    return (R(z + h) - R(z - h)) / (2.0 * h);
}

/*!
 * \brief Trace the stability boundary |R(z)| = 1 by continuation in theta
 *
 * The boundary satisfies R(z(theta)) = exp(i*theta). Starting from the trivial
 * point z(0) = 0 (R(0) = 1 for any consistent method), Newton's method with the
 * previous point as initial guess is used to follow the curve as theta increases.
 *
 * For an explicit Runge-Kutta method of order p, R is a degree-p polynomial whose
 * p roots all lie inside the bounded stability region. By the argument principle,
 * R therefore winds p times around the origin as z traverses the boundary once, so
 * theta has to range over [0, 2*pi*p] for the traced curve to close.
 *
 * The step size in theta is adapted so that consecutive points are roughly
 * `targetStep` apart in the z-plane, taking smaller steps where the curve bends
 * sharply and larger steps along its smoother parts.
 */
std::vector<Complex> traceStabilityBoundary(const ComplexFunction& R, const std::size_t order,
                                             const int stepsPerTurn = 360, const Scalar targetStep = 0.01)
{
    std::vector<Complex> boundary;
    Complex z{0.0, 0.0};
    boundary.push_back(z);

    const Scalar thetaMax = 2.0 * M_PI * order;
    const Scalar dthetaInit = 2.0 * M_PI / stepsPerTurn;
    const Scalar dthetaMin = dthetaInit / 100.0;
    const Scalar dthetaMax = dthetaInit * 4.0;

    Scalar theta = 0.0;
    Scalar dtheta = dthetaInit;
    while (theta < thetaMax - 1e-12)
    {
        const Scalar thetaNext = std::min(theta + dtheta, thetaMax);
        const Complex target = std::polar(Scalar{1.0}, thetaNext);

        Complex zNext = z;
        for (int iter = 0; iter < 50; ++iter)
        {
            const Complex residual = R(zNext) - target;
            if (std::abs(residual) < 1e-13)
                break;
            zNext -= residual / derivative(R, zNext);
        }

        const auto moved = std::abs(zNext - z);
        if (moved > 1e-14)
            dtheta *= std::clamp(targetStep / moved, 0.5, 2.0);
        dtheta = std::clamp(dtheta, dthetaMin, dthetaMax);

        z = zNext;
        theta = thetaNext;
        boundary.push_back(z);
    }

    return boundary;
}

//! Area enclosed by a closed polygon (shoelace formula)
Scalar polygonArea(const std::vector<Complex>& polygon)
{
    Scalar area = 0.0;
    for (std::size_t i = 0; i < polygon.size(); ++i)
    {
        const auto& p0 = polygon[i];
        const auto& p1 = polygon[(i+1) % polygon.size()];
        area += p0.real()*p1.imag() - p1.real()*p0.imag();
    }

    return 0.5 * std::abs(area);
}

//! Largest |Re(z)| at which the traced boundary crosses the negative real axis
Scalar negativeRealAxisCutoff(const std::vector<Complex>& boundary)
{
    Scalar cutoff = 0.0;
    for (std::size_t i = 0; i + 1 < boundary.size(); ++i)
    {
        const auto& p0 = boundary[i];
        const auto& p1 = boundary[i+1];
        if ((p0.imag() >= 0.0) != (p1.imag() >= 0.0))
        {
            const auto t = p0.imag() / (p0.imag() - p1.imag());
            const auto re = p0.real() + t * (p1.real() - p0.real());
            if (re < 0.0)
                cutoff = std::max(cutoff, -re);
        }
    }

    return cutoff;
}

Scalar sampleLeftHalfPlaneMaxModulus(const Method& method,
                                     const Scalar realMin = -6.0,
                                     const Scalar realMax = 0.0,
                                     const Scalar imagAbsMax = 6.0,
                                     const int numReal = 121,
                                     const int numImag = 121)
{
    Scalar maxModulus = 0.0;
    for (int i = 0; i < numReal; ++i)
    {
        const auto real = realMin + (realMax - realMin) * static_cast<Scalar>(i) / (numReal - 1);
        for (int j = 0; j < numImag; ++j)
        {
            const auto imag = -imagAbsMax + 2.0 * imagAbsMax * static_cast<Scalar>(j) / (numImag - 1);
            maxModulus = std::max(maxModulus, std::abs(stabilityFunction(method, Complex{real, imag})));
        }
    }

    return maxModulus;
}

//! Sample points inside the stability region |R(z)| <= 1 on a grid (for plotting unbounded/implicit regions)
std::vector<std::array<Scalar, 2>> sampleStablePoints(const Method& method,
                                                     const Scalar limit = 4.0,
                                                     const int resolution = 121)
{
    std::vector<std::array<Scalar, 2>> points;
    points.reserve(resolution*resolution/2);

    for (int i = 0; i < resolution; ++i)
    {
        const auto imag = -limit + 2.0 * limit * static_cast<Scalar>(i) / (resolution - 1);
        for (int j = 0; j < resolution; ++j)
        {
            const auto real = -limit + 2.0 * limit * static_cast<Scalar>(j) / (resolution - 1);
            if (std::abs(stabilityFunction(method, Complex{real, imag})) <= 1.0 + 1e-12)
                points.push_back({real, imag});
        }
    }

    return points;
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Dumux::Experimental::MultiStage;

    Dumux::initialize(argc, argv);

    const Method explicitEuler = std::make_shared<ExplicitEuler<Scalar>>();
    const Method implicitEuler = std::make_shared<ImplicitEuler<Scalar>>();
    const Method crankNicolson = std::make_shared<Theta<Scalar>>(0.5);
    const Method heun = std::make_shared<RungeKuttaExplicitSecondOrderHeun<Scalar>>();
    const Method rk3 = std::make_shared<RungeKuttaExplicitThirdOrder<Scalar>>();
    const Method rk4 = std::make_shared<RungeKuttaExplicitFourthOrder<Scalar>>();
    const Method dirk2 = std::make_shared<DIRKSecondOrderAlexander<Scalar>>();
    const Method dirk3 = std::make_shared<DIRKThirdOrderAlexander<Scalar>>();

    // shorter display names for plot legends/titles
    const auto displayName = [](const std::string& id, const Method& method)
    {
        if (id == "crank_nicolson") return std::string{"Crank-Nicolson"};
        if (id == "dirk2") return std::string{"DIRK2 (Alexander)"};
        if (id == "dirk3") return std::string{"DIRK3 (Alexander)"};
        return method->name();
    };

    // single-point sanity checks for the stability functions of the implicit methods
    const Complex z(-0.5, 0.0);
    expectNear(std::abs(stabilityFunction(implicitEuler, z) - (1.0 / (1.0 - z))), 0.0, 1e-14,
               "Unexpected stability function for implicit Euler");
    expectNear(std::abs(stabilityFunction(crankNicolson, z) - ((1.0 + 0.5 * z) / (1.0 - 0.5 * z))), 0.0, 1e-14,
               "Unexpected stability function for Crank-Nicolson");

    const Scalar gamma = 1.0 - 1.0 / std::sqrt(2.0);
    expectNear(std::abs(stabilityFunction(dirk2, z)
                        - ((1.0 + (1.0 - 2.0 * gamma) * z) / ((1.0 - gamma * z) * (1.0 - gamma * z)))),
               0.0, 1e-14,
               "Unexpected stability function for Alexander DIRK2");

    // --- explicit methods: trace the stability boundary by continuation and compare
    //     against the analytical truncated-exponential stability polynomial ---
    struct ExplicitMethodData
    {
        std::string id;
        Method method;
        std::size_t order;
        Scalar referenceCutoff;
    };

    const std::vector<ExplicitMethodData> explicitMethods = {
        {"explicit_euler", explicitEuler, 1, 2.0},
        {"heun",           heun,          2, 2.0},
        {"rk3",            rk3,           3, 2.5127},
        {"rk4",            rk4,           4, 2.7853}
    };

    std::vector<std::vector<Complex>> boundaries;
    std::vector<std::vector<Complex>> boundariesAnalytical;
    std::vector<Scalar> areas;

    std::cout << "Negative real-axis cutoffs from traced stability boundaries:" << std::endl;
    for (const auto& entry : explicitMethods)
    {
        const ComplexFunction R = [&](const Complex zz){ return stabilityFunction(entry.method, zz); };
        const ComplexFunction RAnalytical = [&](const Complex zz){ return truncatedExponential(zz, entry.order); };

        auto boundary = traceStabilityBoundary(R, entry.order);
        auto boundaryAnalytical = traceStabilityBoundary(RAnalytical, entry.order);

        expectNear(std::abs(boundary.back() - boundary.front()), 0.0, 1e-6,
                   "Stability boundary continuation did not close for " + entry.method->name());

        for (const auto& point : boundary)
        {
            expectNear(std::abs(R(point)), 1.0, 1e-8,
                       "Traced point is not on the stability boundary |R(z)|=1 for " + entry.method->name());
            expectNear(std::abs(R(point) - truncatedExponential(point, entry.order)), 0.0, 1e-10,
                       "Stability function does not match the analytical truncated exponential for " + entry.method->name());
        }

        const auto cutoff = negativeRealAxisCutoff(boundary);
        std::cout << "  " << entry.method->name() << ": " << cutoff << std::endl;
        expectNear(cutoff, entry.referenceCutoff, 1e-3,
                   "Unexpected negative real-axis cutoff for " + entry.method->name());

        areas.push_back(polygonArea(boundary));
        boundaries.push_back(std::move(boundary));
        boundariesAnalytical.push_back(std::move(boundaryAnalytical));
    }

    // higher-order explicit methods have larger stability regions
    for (std::size_t i = 0; i + 1 < areas.size(); ++i)
        expect(areas[i] < areas[i+1],
               explicitMethods[i+1].method->name() + " should have a larger stability region than "
               + explicitMethods[i].method->name());

    // --- implicit methods: A-stability check in the left half plane ---
    for (const auto& method : {implicitEuler, crankNicolson, dirk2, dirk3})
    {
        const auto maxModulus = sampleLeftHalfPlaneMaxModulus(method);
        expect(maxModulus <= 1.0 + 1e-10,
               "Sampled left-half-plane stability violated for " + method->name());
    }

    const Complex stiffZ(-1.0e6, 0.0);
    expect(std::abs(stabilityFunction(implicitEuler, stiffZ)) < 1.0e-5,
           "Implicit Euler should damp stiff modes");
    expect(std::abs(stabilityFunction(dirk2, stiffZ)) < 1.0e-5,
           "Alexander DIRK2 should damp stiff modes");
    expect(std::abs(stabilityFunction(dirk3, stiffZ)) < 1.0e-5,
           "Alexander DIRK3 should damp stiff modes");
    expect(std::abs(std::abs(stabilityFunction(crankNicolson, stiffZ)) - 1.0) < 1.0e-5,
           "Crank-Nicolson should not be L-stable");

    // --- write data for plotting ---
    const Scalar plotLimit = 4.0;
    const int plotResolution = 121;

    Dumux::Json::JsonTree out;
    out["limit"] = plotLimit;
    out["resolution"] = plotResolution;
    out["methods"] = Dumux::Json::JsonTree::object();

    for (std::size_t i = 0; i < explicitMethods.size(); ++i)
    {
        auto& node = out["methods"][explicitMethods[i].id];
        node["name"] = displayName(explicitMethods[i].id, explicitMethods[i].method);
        node["negativeRealAxisCutoff"] = negativeRealAxisCutoff(boundaries[i]);

        std::vector<std::array<Scalar, 2>> boundaryPoints;
        boundaryPoints.reserve(boundaries[i].size());
        for (const auto& point : boundaries[i])
            boundaryPoints.push_back({point.real(), point.imag()});
        node["boundaryPoints"] = boundaryPoints;

        std::vector<std::array<Scalar, 2>> boundaryPointsAnalytical;
        boundaryPointsAnalytical.reserve(boundariesAnalytical[i].size());
        for (const auto& point : boundariesAnalytical[i])
            boundaryPointsAnalytical.push_back({point.real(), point.imag()});
        node["boundaryPointsAnalytical"] = boundaryPointsAnalytical;
    }

    const std::vector<std::pair<std::string, Method>> implicitMethodsForPlot = {
        {"implicit_euler", implicitEuler},
        {"crank_nicolson", crankNicolson},
        {"dirk2", dirk2},
        {"dirk3", dirk3}
    };

    for (const auto& [id, method] : implicitMethodsForPlot)
    {
        auto& node = out["methods"][id];
        node["name"] = displayName(id, method);
        node["stablePoints"] = sampleStablePoints(method, plotLimit, plotResolution);
    }

    const std::string outputFile = "test_timestepmethods_stabilityregions_data.json";
    std::ofstream output(outputFile);
    output << std::setw(2) << out << std::endl;
    std::cout << "Wrote stability data to " << outputFile << std::endl;

    return 0;
}
