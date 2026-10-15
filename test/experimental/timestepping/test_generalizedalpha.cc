// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Test of the generalized-alpha time stepper for second time derivatives
 *
 * The test system is M u'' + C u' + K u + k3 u^3 = 0 with two degrees of freedom, i.e. the
 * inertia operator M u, the storage operator C u and the spatial operator K u + k3 u^3. It checks
 * - the parameters of the Chung-Hulbert, Hilber-Hughes-Taylor, Wood-Bossak-Zienkiewicz and Newmark methods,
 * - that the stepper (which advances operator values) reproduces the generalized-alpha method
 *   written for the displacement, velocity and acceleration for non-proportional damping and a
 *   nonlinear stiffness,
 * - the observed order for two damped linear oscillators with exact solution,
 * - spectral radius, algorithmic damping ratio and period elongation of the amplification matrix
 *   of the stepper for the undamped oscillator, against the exact values for the Newmark
 *   average acceleration method,
 * - second order for the first-order equation du/dt = -u + exp(t) (storage operator only),
 * - the damping of the constraint of u0' = -u0 + u1, 0 = sin(t) - u1 by -alpha_f/(1 - alpha_f) per step,
 * - that the Newmark parameters reproduce Experimental::NewmarkBeta.
 */
#include <config.h>

#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <memory>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <map>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/io/format.hh>
#include <dumux/io/json.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/experimental/common/variables.hh>
#include <dumux/experimental/ode/odesolver.hh>
#include <dumux/experimental/timestepping/newmarkbeta.hh>
#include <dumux/experimental/timestepping/generalizedalpha.hh>

#include "scalarassembler.hh"
#include "constrainedode.hh"
#include "polynomialroots.hh"

namespace Dumux {

using Scalar = double;
using Vector = Dune::FieldVector<Scalar, 2>;
using Matrix = Dune::FieldMatrix<Scalar, 2, 2>;

void expect(const bool condition, const std::string& message)
{
    if (!condition)
        DUNE_THROW(Dune::InvalidStateException, message);
}

//! M u'' + C u' + K u + k3 u^3 = 0 (cubic term componentwise)
class OscillatorAssembler
{
public:
    using Scalar = double;
    using SolutionVector = Vector;
    using ResidualType = Vector;
    using JacobianMatrix = Matrix;
    using Variables = Experimental::Variables<SolutionVector>;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, ResidualType>;

    OscillatorAssembler(const Matrix& mass, const Matrix& damping, const Matrix& stiffness, const Scalar cubicStiffness)
    : mass_(mass), damping_(damping), stiffness_(stiffness), cubicStiffness_(cubicStiffness)
    {}

    static constexpr bool hasInertia()
    { return true; }

    void setLinearSystem() {}
    JacobianMatrix& jacobian() { return jac_; }
    ResidualType& residual() { return res_; }
    void resetTimeStep(const SolutionVector&) {}

    void setImplicitSolve(Variables& vars, const SolveParams& params)
    {
        vars.updateTime(params.timeLevel);
        params_ = params;
    }

    void assembleOperators(Variables& vars,
                           const Experimental::TimeLevel<Scalar>& timeLevel,
                           Experimental::OperatorEvaluations<ResidualType> evaluations)
    {
        vars.updateTime(timeLevel);
        if (evaluations.inertia)
            mass_.mv(vars.dofs(), *evaluations.inertia);
        if (evaluations.storage)
            damping_.mv(vars.dofs(), *evaluations.storage);
        if (evaluations.spatial)
            *evaluations.spatial = spatial(vars.dofs());
    }

    void assembleResidual(const Variables& vars)
    {
        const auto& w = params_->weights;
        const auto& u = vars.dofs();
        ResidualType inertia, storage;
        mass_.mv(u, inertia);
        damping_.mv(u, storage);
        res_ = 0.0;
        res_.axpy(w.inertia, inertia);
        res_.axpy(w.storage, storage);
        res_.axpy(w.spatial, spatial(u));
        if (params_->constant)
            res_ += *params_->constant;
    }

    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);
        const auto& w = params_->weights;
        jac_ = 0.0;
        jac_.axpy(w.inertia, mass_);
        jac_.axpy(w.storage, damping_);
        jac_.axpy(w.spatial, spatialJacobian(vars.dofs()));
    }

    ResidualType spatial(const Vector& u) const
    {
        ResidualType r;
        stiffness_.mv(u, r);
        for (int i = 0; i < 2; ++i)
            r[i] += cubicStiffness_*u[i]*u[i]*u[i];
        return r;
    }

    Matrix spatialJacobian(const Vector& u) const
    {
        Matrix j = stiffness_;
        for (int i = 0; i < 2; ++i)
            j[i][i] += 3.0*cubicStiffness_*u[i]*u[i];
        return j;
    }

    const Matrix& mass() const { return mass_; }
    const Matrix& damping() const { return damping_; }

private:
    Matrix mass_, damping_, stiffness_;
    Scalar cubicStiffness_;
    ResidualType res_;
    JacobianMatrix jac_;
    std::optional<SolveParams> params_;
};

using LinearSolver = Experimental::ODELocalLinearSolver<Matrix, Vector>;
using Newton = NewtonSolver<OscillatorAssembler, LinearSolver, DefaultPartialReassembler>;
using Method = Experimental::GeneralizedAlphaMethod<Scalar>;

std::shared_ptr<Newton> makeNewton(std::shared_ptr<OscillatorAssembler> assembler)
{
    return std::make_shared<Newton>(assembler, std::make_shared<LinearSolver>(),
                                    Dune::MPIHelper::getCommunication(), "", "Newton", 0);
}

//! The generalized-alpha method for the displacement, velocity and acceleration
struct StateSpaceGeneralizedAlpha
{
    const OscillatorAssembler& system;
    const Method& method;
    Vector u, v, a;

    void step(const Scalar dt)
    {
        const auto am = method.alphaM(), af = method.alphaF(), beta = method.beta(), gamma = method.gamma();
        const auto accelerationAt = [&](const Vector& u1)
        {
            Vector a1 = u1;
            a1 -= u; a1.axpy(-dt, v);
            a1 /= beta*dt*dt;
            a1.axpy(-(0.5 - beta)/beta, a);
            return a1;
        };

        Vector u1 = u;
        for (int iteration = 0; iteration < 50; ++iteration)
        {
            const auto a1 = accelerationAt(u1);
            auto v1 = v; v1.axpy(dt*(1.0 - gamma), a); v1.axpy(dt*gamma, a1);

            Vector accelerationTerm = a1; accelerationTerm *= 1.0 - am; accelerationTerm.axpy(am, a);
            Vector velocityTerm = v1; velocityTerm *= 1.0 - af; velocityTerm.axpy(af, v);
            Vector residual, tmp;
            system.mass().mv(accelerationTerm, residual);
            system.damping().mv(velocityTerm, tmp); residual += tmp;
            residual.axpy(1.0 - af, system.spatial(u1));
            residual.axpy(af, system.spatial(u));

            Matrix jacobian(0.0);
            jacobian.axpy((1.0 - am)/(beta*dt*dt), system.mass());
            jacobian.axpy((1.0 - af)*gamma/(beta*dt), system.damping());
            jacobian.axpy(1.0 - af, system.spatialJacobian(u1));

            Vector update;
            jacobian.solve(update, residual);
            u1 -= update;
            if (update.infinity_norm() < 1e-15)
                break;
        }

        const auto a1 = accelerationAt(u1);
        v.axpy(dt*(1.0 - gamma), a); v.axpy(dt*gamma, a1);
        a = a1;
        u = u1;
    }
};

//! The eigenvalues of the amplification matrix of the stepper for u'' + omega^2 u = 0 with dt = 1
std::vector<std::complex<Scalar>> amplificationEigenvalues(std::shared_ptr<const Method> method, const Scalar omegaDt)
{
    const Matrix mass{{1.0, 0.0}, {0.0, 1.0}};
    const Matrix stiffness{{omegaDt*omegaDt, 0.0}, {0.0, 1.0}};
    auto assembler = std::make_shared<OscillatorAssembler>(mass, Matrix(0.0), stiffness, 0.0);

    // the columns are the states after one step from the unit states (u, v, a)
    Dune::FieldMatrix<Scalar, 3, 3> A;
    for (int j = 0; j < 3; ++j)
    {
        Experimental::GeneralizedAlphaTimeStepper<Newton> stepper(makeNewton(assembler), method);
        OscillatorAssembler::Variables vars(Vector{j == 0 ? 1.0 : 0.0, 0.0});
        stepper.initialize(vars, 0.0, Vector{j == 1 ? 1.0 : 0.0, 0.0}, Vector{j == 2 ? 1.0 : 0.0, 0.0});
        stepper.step(vars, 0.0, 1.0);
        const auto& history = stepper.inertiaHistory();
        A[0][j] = history.value[0];
        A[1][j] = history.rate[0];
        A[2][j] = history.secondRate[0];
    }

    const Scalar trace = A[0][0] + A[1][1] + A[2][2];
    const Scalar minors = A[0][0]*A[1][1] - A[0][1]*A[1][0] + A[0][0]*A[2][2] - A[0][2]*A[2][0]
                        + A[1][1]*A[2][2] - A[1][2]*A[2][1];
    return Test::roots({1.0, -trace, minors, -A.determinant()});
}

//! The eigenvalue r exp(i theta) with theta > 0 that approximates exp(i omega dt)
std::optional<std::complex<Scalar>> principalEigenvalue(const std::vector<std::complex<Scalar>>& eigenvalues)
{
    for (const auto& lambda : eigenvalues)
        if (lambda.imag() > 1e-12)
            return lambda;
    return {};
}

//! The algorithmic damping ratio -ln(r)/theta of the principal eigenvalue r exp(i theta)
Scalar dampingRatio(const std::complex<Scalar>& principal)
{ return -std::log(std::abs(principal))/std::arg(principal); }

//! The ratio of the numerical to the exact period, omega dt/theta
Scalar periodRatio(const std::complex<Scalar>& principal, const Scalar omegaDt)
{ return omegaDt/std::arg(principal); }

//! The exact solution of u'' + c u' + k u = 0 (underdamped)
Scalar dampedOscillator(const Scalar t, const Scalar c, const Scalar k, const Scalar u0, const Scalar v0)
{
    const auto omega = std::sqrt(k);
    const auto zeta = c/(2.0*omega);
    const auto omegaD = omega*std::sqrt(1.0 - zeta*zeta);
    return std::exp(-zeta*omega*t)*(u0*std::cos(omegaD*t) + (v0 + zeta*omega*u0)/omegaD*std::sin(omegaD*t));
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv, [](Dune::ParameterTree& params) {
        params["Newton.MaxRelativeShift"] = "1e-13";
    });

    using namespace Experimental::GeneralizedAlpha;

    ////////////////////////////////////////////////////////////
    // parameters
    ////////////////////////////////////////////////////////////
    {
        const auto near = [](Scalar a, Scalar b) { return std::abs(a - b) < 1e-14; };
        const ChungHulbert<Scalar> ch(0.8);
        expect(near(ch.alphaM(), 0.6/1.8) && near(ch.alphaF(), 0.8/1.8), "Chung-Hulbert alpha_m, alpha_f");
        expect(near(ch.gamma(), 0.5 - ch.alphaM() + ch.alphaF()), "Chung-Hulbert gamma");
        expect(near(ch.beta(), 0.25*(1.0 - ch.alphaM() + ch.alphaF())*(1.0 - ch.alphaM() + ch.alphaF())), "Chung-Hulbert beta");
        const HilberHughesTaylor<Scalar> hht(0.8);
        expect(near(hht.alphaM(), 0.0) && near(hht.alphaF(), 0.2/1.8), "Hilber-Hughes-Taylor alpha_m, alpha_f");
        const WoodBossakZienkiewicz<Scalar> wbz(0.8);
        expect(near(wbz.alphaM(), -0.2/1.8) && near(wbz.alphaF(), 0.0), "Wood-Bossak-Zienkiewicz alpha_m, alpha_f");
        const Newmark<Scalar> newmark;
        expect(near(newmark.beta(), 0.25) && near(newmark.gamma(), 0.5), "Newmark default parameters");

        bool thrown = false;
        try { HilberHughesTaylor<Scalar>(0.4); } catch (const Dune::InvalidStateException&) { thrown = true; }
        expect(thrown, "Hilber-Hughes-Taylor requires rho_inf >= 1/2");
    }

    struct MethodInfo { std::string id; std::shared_ptr<const Method> method; };
    const std::vector<MethodInfo> methods = {
        {"newmark", std::make_shared<Newmark<Scalar>>()},
        {"ch_1.0", std::make_shared<ChungHulbert<Scalar>>(1.0)},
        {"ch_0.8", std::make_shared<ChungHulbert<Scalar>>(0.8)},
        {"ch_0.5", std::make_shared<ChungHulbert<Scalar>>(0.5)},
        {"ch_0.0", std::make_shared<ChungHulbert<Scalar>>(0.0)},
        {"hht_0.8", std::make_shared<HilberHughesTaylor<Scalar>>(0.8)},
        {"wbz_0.8", std::make_shared<WoodBossakZienkiewicz<Scalar>>(0.8)}
    };

    Dumux::Json::JsonTree out;
    out["methods"] = Dumux::Json::JsonTree::object();

    ////////////////////////////////////////////////////////////
    // equality with the method for displacement, velocity and acceleration
    ////////////////////////////////////////////////////////////
    {
        const Matrix mass{{2.0, 0.3}, {0.3, 1.0}};
        const Matrix damping{{0.4, -0.1}, {0.05, 0.2}};
        const Matrix stiffness{{3.0, -1.0}, {-1.0, 2.0}};
        auto assembler = std::make_shared<OscillatorAssembler>(mass, damping, stiffness, 0.5);

        const Vector u0{1.0, -0.5}, v0{0.0, 0.3};
        Vector a0, force = assembler->spatial(u0), dampingForce;
        damping.mv(v0, dampingForce);
        force += dampingForce;
        force *= -1.0;
        mass.solve(a0, force);

        for (const auto& [id, method] : methods)
        {
            Experimental::GeneralizedAlphaTimeStepper<Newton> stepper(makeNewton(assembler), method);
            OscillatorAssembler::Variables vars(u0);
            stepper.initialize(vars, 0.0, v0, a0);

            StateSpaceGeneralizedAlpha reference{*assembler, *method, u0, v0, a0};

            const Scalar dt = 0.05;
            Scalar maxDifference = 0.0;
            for (int n = 0; n < 60; ++n)
            {
                stepper.step(vars, n*dt, dt);
                reference.step(dt);
                auto difference = vars.dofs(); difference -= reference.u;
                maxDifference = std::max(maxDifference, difference.infinity_norm());
            }
            std::cout << Fmt::format("-- {}: largest difference to the state-space method {:.2e}\n", id, maxDifference);
            expect(maxDifference < 1e-11, "The stepper does not reproduce the state-space generalized-alpha method for " + id);
        }
    }

    ////////////////////////////////////////////////////////////
    // observed order for damped linear oscillators
    ////////////////////////////////////////////////////////////
    {
        const Scalar c0 = 0.2, k0 = 4.0, c1 = 0.5, k1 = 9.0;
        const Matrix mass{{1.0, 0.0}, {0.0, 1.0}};
        const Matrix damping{{c0, 0.0}, {0.0, c1}};
        const Matrix stiffness{{k0, 0.0}, {0.0, k1}};
        auto assembler = std::make_shared<OscillatorAssembler>(mass, damping, stiffness, 0.0);
        const Vector u0{1.0, 0.5}, v0{0.0, 1.0};
        const Vector a0{-c0*v0[0] - k0*u0[0], -c1*v0[1] - k1*u0[1]};
        const Scalar tEnd = 2.0;

        for (const auto& [id, method] : methods)
        {
            std::vector<std::array<Scalar, 2>> convergence;
            for (const int n : {25, 50, 100, 200, 400})
            {
                Experimental::GeneralizedAlphaTimeStepper<Newton> stepper(makeNewton(assembler), method);
                OscillatorAssembler::Variables vars(u0);
                stepper.initialize(vars, 0.0, v0, a0);
                const Scalar dt = tEnd/n;
                for (int i = 0; i < n; ++i)
                    stepper.step(vars, i*dt, dt);
                const Scalar error = std::max(std::abs(vars.dofs()[0] - dampedOscillator(tEnd, c0, k0, u0[0], v0[0])),
                                              std::abs(vars.dofs()[1] - dampedOscillator(tEnd, c1, k1, u0[1], v0[1])));
                convergence.push_back({dt, error});
            }
            const auto rate = std::log2(convergence[3][1]/convergence[4][1]);
            std::cout << Fmt::format("-- {}: observed order {:.3f}\n", id, rate);
            expect(rate > 1.9, "Observed order of " + id);

            auto& node = out["methods"][id];
            node["name"] = method->name();
            node["convergence"] = convergence;
        }
    }

    ////////////////////////////////////////////////////////////
    // amplification matrix for the undamped oscillator
    ////////////////////////////////////////////////////////////
    for (const auto& [id, method] : methods)
    {
        std::vector<std::array<Scalar, 2>> spectralRadii;
        std::vector<std::array<Scalar, 3>> dampingAndPeriod;
        Scalar maxSpectralRadius = 0.0;
        Scalar stiffSpectralRadius = 0.0;
        for (int i = 0; i <= 240; ++i)
        {
            const Scalar omegaDt = std::pow(10.0, -2.0 + 8.0*i/240.0);
            const auto eigenvalues = amplificationEigenvalues(method, omegaDt);

            Scalar spectralRadius = 0.0;
            for (const auto& lambda : eigenvalues)
                spectralRadius = std::max(spectralRadius, std::abs(lambda));
            maxSpectralRadius = std::max(maxSpectralRadius, spectralRadius);
            stiffSpectralRadius = spectralRadius;

            spectralRadii.push_back({omegaDt, spectralRadius});
            if (i == 0)
                expect(spectralRadius > 0.999, "Low frequencies are not damped by " + id);

            if (const auto principal = principalEigenvalue(eigenvalues))
            {
                const auto damping = dampingRatio(*principal);
                const auto period = periodRatio(*principal, omegaDt);
                dampingAndPeriod.push_back({omegaDt, damping, period});

                if (id == "newmark")
                {
                    expect(std::abs(spectralRadius - 1.0) < 1e-10, "Newmark is not dissipative");
                    expect(std::abs(damping) < 1e-10, "Newmark has no algorithmic damping");
                    const Scalar exactPeriodRatio = omegaDt/(2.0*std::atan(0.5*omegaDt));
                    expect(std::abs(period - exactPeriodRatio) < 1e-8*exactPeriodRatio, "Period of the Newmark method");
                }
            }
        }

        const auto principal = principalEigenvalue(amplificationEigenvalues(method, 0.2*M_PI));
        expect(principal.has_value(), "Complex principal eigenvalue at dt/T = 0.1 of " + id);
        std::cout << Fmt::format("-- {}: at dt/T = 0.1, damping ratio {:.5f}, relative period error {:.5f}\n",
                                 id, dampingRatio(*principal), periodRatio(*principal, 0.2*M_PI) - 1.0);

        // for rho_inf = 1, all three eigenvalues approach -1 for omega dt -> infinity, and the roots
        // of the characteristic polynomial are then only accurate to the cube root of the round-off
        std::cout << Fmt::format("-- {}: largest spectral radius {:.8f}\n", id, maxSpectralRadius);
        expect(maxSpectralRadius < 1.0 + 1e-4, "Unconditional stability of " + id);
        std::cout << Fmt::format("-- {}: spectral radius at omega dt = 1e6: {:.5f}\n", id, stiffSpectralRadius);
        const std::map<std::string, Scalar> expectedRhoInf{
            {"newmark", 1.0}, {"ch_1.0", 1.0}, {"ch_0.8", 0.8}, {"ch_0.5", 0.5}, {"ch_0.0", 0.0}, {"hht_0.8", 0.8}, {"wbz_0.8", 0.8}
        };
        expect(std::abs(stiffSpectralRadius - expectedRhoInf.at(id)) < 1e-3, "Spectral radius at infinity of " + id);

        out["methods"][id]["spectralRadius"] = spectralRadii;
        out["methods"][id]["dampingAndPeriod"] = dampingAndPeriod;
    }

    ////////////////////////////////////////////////////////////
    // first-order equation du/dt = -u + exp(t), u(0) = 1 (storage only)
    ////////////////////////////////////////////////////////////
    for (const auto& [id, method] : methods)
    {
        using ScalarNewton = NewtonSolver<ScalarAssembler, ScalarLinearSolver, DefaultPartialReassembler>;
        std::vector<Scalar> errors;
        for (const int n : {50, 100, 200})
        {
            auto newton = std::make_shared<ScalarNewton>(std::make_shared<ScalarAssembler>(-1.0), std::make_shared<ScalarLinearSolver>(),
                                                         Dune::MPIHelper::getCommunication(), "", "Newton", 0);
            Experimental::GeneralizedAlphaTimeStepper<ScalarNewton> stepper(newton, method);
            ScalarAssembler::Variables vars(1.0);
            stepper.initialize(vars, 0.0, 0.0, 1.0);
            const Scalar dt = 1.0/n;
            for (int i = 0; i < n; ++i)
                stepper.step(vars, i*dt, dt);
            errors.push_back(std::abs(vars.dofs() - std::cosh(1.0)));
        }
        const auto rate = std::log2(errors[1]/errors[2]);
        std::cout << Fmt::format("-- {}: observed order for the first-order equation {:.3f}\n", id, rate);
        expect(rate > 1.9, "Observed order for the first-order equation with " + id);
    }

    ////////////////////////////////////////////////////////////
    // the constraint error of u0' = -u0 + u1, 0 = sin(t) - u1 is multiplied by -alpha_f/(1 - alpha_f)
    ////////////////////////////////////////////////////////////
    for (const auto& [id, method] : methods)
    {
        using ODEAssembler = Experimental::MultiStageODEAssembler<Test::ConstrainedODE>;
        using ODELinearSolver = Experimental::ODELocalLinearSolver<Test::ConstrainedODE::JacobianMatrix, Test::ConstrainedODE::ResidualType>;
        using ODENewton = NewtonSolver<ODEAssembler, ODELinearSolver, DefaultPartialReassembler>;
        auto newton = std::make_shared<ODENewton>(std::make_shared<ODEAssembler>(std::make_shared<const Test::ConstrainedODE>()),
                                                  std::make_shared<ODELinearSolver>(),
                                                  Dune::MPIHelper::getCommunication(), "", "Newton", 0);
        Experimental::GeneralizedAlphaTimeStepper<ODENewton> stepper(newton, method);
        Test::ConstrainedODE::Variables vars(Test::ConstrainedODE::SolutionVector{0.0, 0.1});
        stepper.initialize(vars, 0.0);

        const Scalar factor = -method->alphaF()/(1.0 - method->alphaF());
        Scalar error = 0.1;
        const Scalar dt = 0.1;
        for (int i = 0; i < 4; ++i)
        {
            stepper.step(vars, i*dt, dt);
            const Scalar newError = vars.dofs()[1] - std::sin((i+1)*dt);
            expect(std::abs(newError - factor*error) < 1e-12, "Damping of the constraint error by " + id);
            error = newError;
        }
        std::cout << Fmt::format("-- {}: constraint error after four steps {:.3e} (factor {:.4f})\n", id, error, factor);
    }

    ////////////////////////////////////////////////////////////
    // the Newmark parameters reproduce Experimental::NewmarkBeta for u'' = -u
    ////////////////////////////////////////////////////////////
    {
        const Scalar beta = 0.3, gamma = 0.6, dt = 0.1;
        using Solution = Dune::BlockVector<Scalar>;
        Experimental::NewmarkBeta<Scalar, Solution> newmarkBeta(beta, gamma);
        Solution x(1), v(1), a(1);
        x = -1.0; v = 0.0; a = 1.0;
        newmarkBeta.initialize(x, v, a);

        const Matrix mass{{1.0, 0.0}, {0.0, 1.0}};
        const Matrix stiffness{{1.0, 0.0}, {0.0, 1.0}};
        auto assembler = std::make_shared<OscillatorAssembler>(mass, Matrix(0.0), stiffness, 0.0);
        Experimental::GeneralizedAlphaTimeStepper<Newton> stepper(makeNewton(assembler),
            std::make_shared<Newmark<Scalar>>(beta, gamma));
        OscillatorAssembler::Variables vars(Vector{-1.0, 0.0});
        stepper.initialize(vars, 0.0, Vector{0.0, 0.0}, Vector{1.0, 0.0});

        Scalar maxDifference = 0.0;
        for (int n = 0; n < 100; ++n)
        {
            // the acceleration of NewmarkBeta is linear in the new position: solve a(x) + x = 0
            const Scalar r0 = newmarkBeta.acceleration(0, dt, Solution::block_type(0.0)) + 0.0;
            const Scalar r1 = newmarkBeta.acceleration(0, dt, Solution::block_type(1.0)) + 1.0;
            Solution xNew(1);
            xNew = -r0/(r1 - r0);
            newmarkBeta.update(dt, xNew);

            stepper.step(vars, n*dt, dt);
            maxDifference = std::max(maxDifference, std::abs(vars.dofs()[0] - xNew[0]));
        }
        std::cout << Fmt::format("-- Newmark (beta = {}, gamma = {}): largest difference to NewmarkBeta {:.2e}\n", beta, gamma, maxDifference);
        expect(maxDifference < 1e-12, "The Newmark parameters do not reproduce NewmarkBeta");
    }

    const std::string outputFile = "test_generalizedalpha_data.json";
    std::ofstream output(outputFile);
    output << std::setw(2) << out << std::endl;
    std::cout << "Wrote data to " << outputFile << std::endl;

    return 0;
}
