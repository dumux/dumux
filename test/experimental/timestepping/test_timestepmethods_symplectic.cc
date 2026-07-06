//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Special test for the symplectic Qin-Zhang DIRK scheme (@cite QinZhang1992).
//
// A Runge-Kutta method is symplectic iff its coefficients satisfy the Cooper condition
//   b_i a_ij + b_j a_ji - b_i b_j = 0   for all i, j,
// which is also the condition under which the method conserves every quadratic invariant
// exactly. This test verifies that condition algebraically for the Qin-Zhang tableau and
// then demonstrates its two hallmark consequences numerically on the harmonic oscillator
// (a linear Hamiltonian system):
//   1. the one-step map preserves phase-space area exactly (its determinant is 1), and
//   2. the quadratic energy H = (q^2 + p^2)/2 is conserved to machine precision over many
//      periods, whereas the non-symplectic implicit Euler (1st order) and Alexander DIRK2
//      (same 2nd order) schemes dissipate it.
//
#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/io/format.hh>
#include <dumux/io/json.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>

namespace Dumux {

using Scalar = double;
using Vec2 = std::array<Scalar, 2>;
using Mat2 = std::array<std::array<Scalar, 2>, 2>; // row-major
using Method = std::shared_ptr<const Experimental::MultiStageMethod<Scalar>>;

void expect(const bool cond, const std::string& msg)
{
    if (!cond)
        DUNE_THROW(Dune::InvalidStateException, msg);
}

void expectNear(const Scalar value, const Scalar reference, const Scalar tolerance, const std::string& msg)
{
    using std::abs;
    if (abs(value - reference) > tolerance)
        DUNE_THROW(Dune::InvalidStateException,
                   msg + " (|" + std::to_string(value) + " - " + std::to_string(reference)
                   + "| = " + std::to_string(abs(value - reference)) + " > " + std::to_string(tolerance) + ")");
}

Vec2 matVec(const Mat2& A, const Vec2& x)
{ return {{ A[0][0]*x[0] + A[0][1]*x[1], A[1][0]*x[0] + A[1][1]*x[1] }}; }

//! solve the 2x2 system A y = rhs via Cramer's rule (A is well-conditioned here)
Vec2 solve2(const Mat2& A, const Vec2& rhs)
{
    const Scalar det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    return {{ (rhs[0]*A[1][1] - A[0][1]*rhs[1]) / det,
              (A[0][0]*rhs[1] - rhs[0]*A[1][0]) / det }};
}

Scalar energy(const Vec2& u)
{ return 0.5*(u[0]*u[0] + u[1]*u[1]); }

/*!
 * \brief Advance u' = L u by one step of a multi-stage method, applying its Shu-Osher
 *        recurrence directly to the 2x2 linear system. Being linear, every stage is a
 *        single exact 2x2 solve (no Newton iteration), so the only error is from the
 *        method itself, not from an inexact stage solve.
 *
 * For stage i the residual sum_{k<=i} [alpha_ik M(x^k) + beta_ik dt R(x^k)] = 0 with
 * M(x) = x and R(x) = -L x becomes
 *   (alpha_ii I - dt beta_ii L) x^i = - sum_{k<i} (alpha_ik I - dt beta_ik L) x^k,
 * which is the vector analogue of the scalar stability-function recurrence.
 */
Vec2 advanceLinear(const Method& method, const Scalar dt, const Mat2& L, const Vec2& u0)
{
    const auto m = method->numStages();
    std::vector<Vec2> stage(m + 1, Vec2{{0.0, 0.0}});
    stage[0] = u0;

    // (a I - dt b L) as a 2x2 matrix
    const auto opMatrix = [&](const Scalar a, const Scalar b) -> Mat2 {
        return {{ {{ a - dt*b*L[0][0], -dt*b*L[0][1] }},
                  {{ -dt*b*L[1][0], a - dt*b*L[1][1] }} }};
    };

    for (std::size_t i = 1; i <= m; ++i)
    {
        Vec2 rhs{{0.0, 0.0}};
        for (std::size_t k = 0; k < i; ++k)
        {
            const auto contribution = matVec(opMatrix(method->temporalWeight(i, k),
                                                      method->spatialWeight(i, k)), stage[k]);
            rhs[0] -= contribution[0];
            rhs[1] -= contribution[1];
        }
        stage[i] = solve2(opMatrix(method->temporalWeight(i, i), method->spatialWeight(i, i)), rhs);
    }

    return stage[m];
}

//! determinant of the one-step propagation matrix of u' = L u (columns are the images of e1, e2)
Scalar oneStepDeterminant(const Method& method, const Scalar dt, const Mat2& L)
{
    const auto c0 = advanceLinear(method, dt, L, Vec2{{1.0, 0.0}});
    const auto c1 = advanceLinear(method, dt, L, Vec2{{0.0, 1.0}});
    return c0[0]*c1[1] - c1[0]*c0[1];
}

//! integrate the harmonic oscillator and return {max relative energy drift, final energy / initial energy}
std::pair<Scalar, Scalar> integrateOscillator(const Method& method, const Scalar dt,
                                              const std::size_t numSteps, const Mat2& L)
{
    Vec2 u{{1.0, 0.0}};
    const Scalar H0 = energy(u);

    using std::abs;
    Scalar maxDrift = 0.0;
    for (std::size_t n = 0; n < numSteps; ++n)
    {
        u = advanceLinear(method, dt, L, u);
        maxDrift = std::max(maxDrift, abs(energy(u) - H0)/H0);
    }

    return { maxDrift, energy(u)/H0 };
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Dumux::Experimental::MultiStage;

    Dumux::initialize(argc, argv);

    const Method qinZhang = std::make_shared<QinZhangSymplecticDIRK<Scalar>>();
    const Method implicitEuler = std::make_shared<ImplicitEuler<Scalar>>();
    const Method dirk2 = std::make_shared<DIRKSecondOrderAlexander<Scalar>>();

    // --- 1. extract the Butcher tableau of the Qin-Zhang scheme from its Shu-Osher form ---
    // The scheme is a two-stage DIRK with the (non stiffly-accurate) update x^{n+1} appended
    // as an explicit third stage, so a_ij = beta(i, j), b_j = beta(numStages(), j), c_i = d_i.
    constexpr std::size_t s = 2;
    Mat2 A{{{{0.0, 0.0}}, {{0.0, 0.0}}}};
    Vec2 b{{0.0, 0.0}};
    Vec2 c{{0.0, 0.0}};
    for (std::size_t i = 1; i <= s; ++i)
    {
        c[i-1] = qinZhang->timeStepWeight(i);
        b[i-1] = qinZhang->spatialWeight(qinZhang->numStages(), i);
        for (std::size_t j = 1; j <= s; ++j)
            A[i-1][j-1] = qinZhang->spatialWeight(i, j);
    }

    // the extracted tableau has to match the documented Qin-Zhang coefficients
    constexpr Scalar eps = 1e-14;
    expectNear(A[0][0], 0.25, eps, "Qin-Zhang a_11 should be 1/4");
    expectNear(A[0][1], 0.0,  eps, "Qin-Zhang a_12 should be 0");
    expectNear(A[1][0], 0.5,  eps, "Qin-Zhang a_21 should be 1/2");
    expectNear(A[1][1], 0.25, eps, "Qin-Zhang a_22 should be 1/4");
    expectNear(b[0], 0.5, eps, "Qin-Zhang b_1 should be 1/2");
    expectNear(b[1], 0.5, eps, "Qin-Zhang b_2 should be 1/2");
    expectNear(c[0], 0.25, eps, "Qin-Zhang c_1 should be 1/4");
    expectNear(c[1], 0.75, eps, "Qin-Zhang c_2 should be 3/4");

    // consistency: sum of weights is 1 and the row sums of A equal the abscissae c
    expectNear(b[0] + b[1], 1.0, eps, "Qin-Zhang weights b have to sum to 1");
    for (std::size_t i = 0; i < s; ++i)
        expectNear(A[i][0] + A[i][1], c[i], eps, "Qin-Zhang row sum of A must equal c");

    // --- 2. algebraic symplecticity: Cooper condition b_i a_ij + b_j a_ji - b_i b_j = 0 ---
    std::cout << "-- Symplecticity matrix M_ij = b_i a_ij + b_j a_ji - b_i b_j:\n";
    for (std::size_t i = 0; i < s; ++i)
    {
        for (std::size_t j = 0; j < s; ++j)
        {
            const Scalar mij = b[i]*A[i][j] + b[j]*A[j][i] - b[i]*b[j];
            std::cout << Fmt::format("   M_{}{} = {:+.3e}\n", i+1, j+1, mij);
            expectNear(mij, 0.0, eps,
                       "Qin-Zhang violates the symplecticity condition at (" + std::to_string(i+1)
                       + "," + std::to_string(j+1) + ")");
        }
    }
    std::cout << "-- Qin-Zhang satisfies the symplecticity (Cooper) condition.\n\n";

    // --- 3. numerical benchmark: harmonic oscillator q' = p, p' = -q ---
    // system matrix L for u = (q, p), i.e. u' = L u, with quadratic invariant H = (q^2+p^2)/2
    const Mat2 L{{{{0.0, 1.0}}, {{-1.0, 0.0}}}};

    // 3a. phase-space area preservation: the one-step map determinant is 1 for a symplectic
    //     method (area preserving) but < 1 for the dissipative reference schemes
    const Scalar dtMap = 0.5;
    const auto detQinZhang = oneStepDeterminant(qinZhang, dtMap, L);
    const auto detImplicitEuler = oneStepDeterminant(implicitEuler, dtMap, L);
    const auto detDirk2 = oneStepDeterminant(dirk2, dtMap, L);

    std::cout << "-- One-step map determinant (phase-space area factor) at dt = " << dtMap << ":\n"
              << Fmt::format("   Qin-Zhang:      det = {:.15f}\n", detQinZhang)
              << Fmt::format("   implicit Euler: det = {:.15f}\n", detImplicitEuler)
              << Fmt::format("   DIRK2:          det = {:.15f}\n\n", detDirk2);

    expectNear(detQinZhang, 1.0, 1e-12,
               "Symplectic Qin-Zhang must preserve phase-space area (one-step determinant = 1)");
    expect(detImplicitEuler < 1.0 - 1e-3,
           "implicit Euler should contract phase-space area (determinant < 1)");
    expect(detDirk2 < 1.0 - 1e-5,
           "Alexander DIRK2 should contract phase-space area (determinant < 1)");

    // 3b. long-time energy conservation over 200 periods
    const std::size_t stepsPerPeriod = 50;
    const std::size_t numPeriods = 200;
    const Scalar dtLong = 2.0*M_PI/static_cast<Scalar>(stepsPerPeriod);
    const std::size_t numStepsLong = stepsPerPeriod*numPeriods;

    const auto [driftQinZhang, finalQinZhang] = integrateOscillator(qinZhang, dtLong, numStepsLong, L);
    const auto [driftImplicitEuler, finalImplicitEuler] = integrateOscillator(implicitEuler, dtLong, numStepsLong, L);
    const auto [driftDirk2, finalDirk2] = integrateOscillator(dirk2, dtLong, numStepsLong, L);

    std::cout << "-- Energy over " << numPeriods << " periods (" << numStepsLong << " steps), H_0 = 0.5:\n"
              << Fmt::format("   Qin-Zhang:      max rel. drift = {:.3e}, final E/E_0 = {:.6f}\n", driftQinZhang, finalQinZhang)
              << Fmt::format("   implicit Euler: max rel. drift = {:.3e}, final E/E_0 = {:.6e}\n", driftImplicitEuler, finalImplicitEuler)
              << Fmt::format("   DIRK2:          max rel. drift = {:.3e}, final E/E_0 = {:.6e}\n\n", driftDirk2, finalDirk2);

    // the symplectic scheme conserves the quadratic energy to machine precision ...
    expect(driftQinZhang < 1e-10,
           "Symplectic Qin-Zhang should conserve the quadratic energy to machine precision, but drift = "
           + std::to_string(driftQinZhang));
    // ... while the non-symplectic schemes visibly dissipate it (implicit Euler strongly, the
    //     same-order DIRK2 more slowly) -- so it is symplecticity, not the order, that conserves energy
    expect(finalImplicitEuler < 0.5,
           "implicit Euler should strongly dissipate the oscillator energy");
    expect(driftDirk2 > 1e6*driftQinZhang,
           "the non-symplectic DIRK2 should dissipate energy far more than the symplectic Qin-Zhang");

    // --- 4. write phase-space trajectories and energy histories for plotting ---
    const std::size_t stepsPerPeriodPlot = 60;
    const std::size_t numPeriodsPlot = 8;
    const Scalar dtPlot = 2.0*M_PI/static_cast<Scalar>(stepsPerPeriodPlot);
    const std::size_t numStepsPlot = stepsPerPeriodPlot*numPeriodsPlot;

    const std::vector<std::pair<std::string, Method>> plotMethods = {
        {"qin_zhang", qinZhang},
        {"dirk2", dirk2},
        {"implicit_euler", implicitEuler}
    };
    const auto displayName = [](const std::string& id) {
        if (id == "qin_zhang") return std::string{"Qin-Zhang (symplectic)"};
        if (id == "dirk2") return std::string{"DIRK2 (Alexander)"};
        if (id == "implicit_euler") return std::string{"implicit Euler"};
        return id;
    };

    Dumux::Json::JsonTree out;
    out["dt"] = dtPlot;
    out["periods"] = numPeriodsPlot;
    out["methods"] = Dumux::Json::JsonTree::object();
    for (const auto& [id, method] : plotMethods)
    {
        auto& node = out["methods"][id];
        node["name"] = displayName(id);

        Vec2 u{{1.0, 0.0}};
        const Scalar H0 = energy(u);
        std::vector<std::array<Scalar, 2>> trajectory;
        std::vector<std::array<Scalar, 2>> energyHistory;
        trajectory.reserve(numStepsPlot + 1);
        energyHistory.reserve(numStepsPlot + 1);
        trajectory.push_back({u[0], u[1]});
        energyHistory.push_back({0.0, energy(u)/H0});
        for (std::size_t n = 0; n < numStepsPlot; ++n)
        {
            u = advanceLinear(method, dtPlot, L, u);
            trajectory.push_back({u[0], u[1]});
            energyHistory.push_back({(n + 1)*dtPlot, energy(u)/H0});
        }
        node["trajectory"] = trajectory;
        node["energy"] = energyHistory;
    }

    const std::string outputFile = "test_timestepmethods_symplectic_data.json";
    std::ofstream output(outputFile);
    output << std::setw(2) << out << std::endl;
    std::cout << "Wrote symplecticity benchmark data to " << outputFile << std::endl;

    return 0;
}
