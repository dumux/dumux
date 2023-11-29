// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief MLMC diffusion equation with random scalar diffusion coefficient
 *
 * Test based on https://bitbucket.org/pefarrell/pymlmc/src/master/pymlmc/mlmc_test.py
 * published under GNU GPL (© Farrell, Giles, Croci, Roy, Beentjes).
 */
#include <config.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

#include <dumux/uq/mlmc.hh>

#include "diffusion_solver.hh"

template <class Scalar = double>
class DiffusionScalarProblem
{
public:
    DiffusionScalarProblem(int maxLevel)
    {
        for (int l = 0; l < maxLevel+1; ++l)
            solvers_.emplace_back(l);
    }

    Scalar evaluate(int level, const Scalar D) const
    {
        auto* orig_buf = std::cout.rdbuf();
        std::cout.rdbuf(NULL);

        const auto QoI = solvers_[level].solve(D, 0.1, 30);

        std::cout.rdbuf(orig_buf);

        return QoI;
    }
private:
    std::vector<const Dumux::DiffusionSolver<2>> solvers_;
};

template<class Scalar>
auto sample(int N, int l, std::default_random_engine& gen, std::normal_distribution<Scalar>& dis)
{
    if (N != 1)
        DUNE_THROW(Dune::NotImplemented, "Vectorized sampling");

    using std::exp;
    const auto s = exp(dis(gen)); // log-normal distribution
    return std::make_tuple(s, s);
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv, [](auto& p){
        p["Problem.Name"] = "diff";
        p["LinearSolver.MaxIterations"] = "10000";
    });

    std::default_random_engine gen;
    const double sigma = Dumux::getParam<double>("Sampler.Variance", 1.0);
    std::normal_distribution<double> dis(0.0, std::sqrt(sigma));
    const auto sampler = [&](int N, int l){
        return sample(N, l, gen, dis);
    };

    const auto N0 = getParam<int>("MLMC.InitialNumSamplesOnCoarseLevels", 10);
    const auto Lmin = getParam<int>("MLMC.MinRefinementLevel", 2);
    const auto Lmax = getParam<int>("MLMC.MaxRefinementLevel", 4);
    const auto N = getParam<int>("MLMC.SamplesForConvergenceTest", 100);
    const auto L = getParam<int>("MLMC.LevelsForConvergenceTest", 4);

    // user defined functional
    DiffusionScalarProblem problem(Lmax);
    const auto func = [&](int level, int numSamples) {
        return Dumux::mlmcFunction(level, numSamples, problem, sampler);
    };

    // run tests
    std::cout << "\n"
              << "**********************************************************\n"
              << "*** Convergence tests, kurtosis, telescoping sum check ***\n"
              << "*** using N = " << N << " samples\n"
              << "**********************************************************\n\n"
              << " l  ave(Pf-Pc)   ave(Pf)     var(Pf-Pc)  var(Pf)      kurtosis     check      cost/sample\n"
              << "----------------------------------------------------------------------------------\n";

    std::vector<double> cost, del1, del2, var1, var2, kur1, chk1;
    for (int l = 0; l < L+1; ++l)
    {
        using std::max, std::abs, std::sqrt;

        auto [sums, cst] = func(l, N);
        for (auto& s : sums)
            s /= N;
        cst /= N;

        const auto kurtosis = (l == 0) ? 0.0 : (
            sums[3]
            - 4*sums[2]*sums[0]
            + 6*sums[1]*sums[0]*sums[0]
            - 3*sums[0]*sums[0]*sums[0]*sums[0]
        ) / ((sums[1]-sums[0]*sums[0])*(sums[1]-sums[0]*sums[0]));

        cost.push_back(cst);
        del1.push_back(sums[0]);
        del2.push_back(sums[4]);
        var1.push_back(max(sums[1]-sums[0]*sums[0], 1e-10));
        var2.push_back(max(sums[5]-sums[4]*sums[4], 1e-10)); // fix for cases with var = 0
        kur1.push_back(kurtosis);

        const auto check = [&]{
            if (l == 0) return 0.0;
            return sqrt(double(N)) * abs(del1[l] + del2[l-1] - del2[l]) / (
                3.0*( sqrt(var1[l]) + sqrt(var2[l-1]) + sqrt(var2[l]) )
            );
        }();

        chk1.push_back(check);

        std::cout << Fmt::format(
            "{:2d} {:11.4e} {:11.4e} {:11.3e} {:11.3e} {:11.2e} {:11.2e} {:11.2e}\n",
            l, del1[l], del2[l], var1[l], var2[l], kur1[l], chk1[l], cst
        );
    }

    if (kur1.back() > 100.0)
        std::cout << Fmt::format("\n WARNING: kurtosis on finest level = {} \n", kur1[-1])
                  << " indicates MLMC correction dominated by a few rare paths; \n"
                  << " for information on the connection to variance of sample variances,\n"
                  << " see http://mathworld.wolfram.com/SampleVarianceDistribution.html\n\n";

    if (auto max = *std::max_element(chk1.begin(), chk1.end()); max > 1.0)
        std::cout << Fmt::format("\n WARNING: maximum consistency error = {} \n", max)
                  << " indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied; \n"
                  << " to be more certain, re-run mlmc_test with larger N \n\n";

    // use linear regression to estimate alpha, beta, gamma
    std::vector<double> x(L, 0.0), y(L, 0.0);
    for (int l=1; l<=L; ++l)
    {
        x[l-1] = l;
        y[l-1] = -std::log2(std::abs(del1[l]));
    }
    const auto [i0, alpha] = linearRegression(x, y);

    for (int l=1; l<=L; ++l)
    {
        x[l-1] = l;
        y[l-1] = -std::log2(std::abs(var1[l]));
    }
    const auto [i1, beta] = linearRegression(x, y);

    for (int l=1; l<=L; ++l)
    {
        x[l-1] = l;
        y[l-1] = -std::log2(std::abs(cost[l]));
    }
    const auto [i2, gamma] = linearRegression(x, y);

    std::cout << "\n"
              << "******************************************************\n"
              << "*** Linear regression estimates of MLMC parameters ***\n"
              << "******************************************************\n"
              << Fmt::format("\n alpha = {:6.3f} (exponent for MLMC weak convergence)\n", alpha)
              << Fmt::format(" beta  = {:6.3f} (exponent for MLMC variance) \n", beta)
              << Fmt::format(" gamma = {:6.3f} (exponent for MLMC cost) \n", -gamma);


    std::cout << "\n"
              << "*****************************\n"
              << "*** MLMC complexity tests ***\n"
              << "*****************************\n\n"
              << "   eps      value     MLMC cost   MC cost  savings     N_l\n"
              << "---------------------------------------------------------------------------------\n";

    std::vector<double> epsilons{ 0.002, 0.005, 0.01, 0.02, 0.05, 0.1 };

    using std::max, std::min;
    const auto alpha1 = max(alpha, 0.5);
    const auto beta1  = max(beta, 0.5);
    const auto gamma1  = -gamma;
    const auto theta = 0.25;

    for (auto eps : epsilons)
    {
        const auto [P, Nl, Cl] = mlmc(Lmin, Lmax, N0, eps, func, alpha1, beta1, gamma1);
        const auto maxLevel = Nl.size() - 1;

        const auto mlmcCost = std::inner_product(Nl.begin(), Nl.end(), Cl.begin(), 1.0);
        const auto mcCost  = var2.back()*Cl[min(Cl.size()-1, maxLevel)]/((1.0-theta)*eps*eps);

        std::cout << Fmt::format("{:5.3e}  {:11.4e}  {:5.3e}  {:5.3e}  {:7.2f}", eps, P, mlmcCost, mcCost, mcCost/mlmcCost);
        for (int l = 0; l < Nl.size(); ++l)
            std::cout << Fmt::format("{:8d}", Nl[l]);
        std::cout << "\n";
    }

    return 0;
}
