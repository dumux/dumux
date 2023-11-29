// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup UncertaintyQuantification
 * \brief Multi-level Monte Carlo sampler
 *
 * Based on details in Giles 2014 (10.1017/S096249291500001X)
 * and the Python implementation of Farrell, Giles, Croci, Roy, Beentjes
 * (https://bitbucket.org/pefarrell/pymlmc) released under GNU GPL.
 */
#ifndef DUMUX_UQ_MLMC_HH
#define DUMUX_UQ_MLMC_HH

#include <chrono>

#include <dumux/common/exceptions.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup UncertaintyQuantification
 * \brief MLMC functional
 *
 * TODO this choice of paramters seems a bit odd. We only
 * need the l and l-1 problems (or only problem 0 for l==0).
 * Number of problems can be deduced from the problem type.
 *
 * Inputs:
        l: level
        N: number of paths
        problems: list of problems
            problems[l-1]: application-specific coarse problem (for l>0)
            problems[l]: application-specific fine problem
            Problems must have an evaluate method such that
            problems[l].evaluate(sample) returns output P_l.
            Optionally, user-defined problems.cost
        sampler: sampling function, by default standard Normal.
            input: N, l
            output: (samplef, samplec). The fine and coarse samples.

    Outputs:
        (sums, cost) where sums is an array of outputs:
        sums[0] = sum(Pf-Pc)
        sums[1] = sum((Pf-Pc)**2)
        sums[2] = sum((Pf-Pc)**3)
        sums[3] = sum((Pf-Pc)**4)
        sums[4] = sum(Pf)
        sums[5] = sum(Pf**2)
        cost = user-defined computational cost. By default, time
 */
template<class Problem, class Sampler, class Scalar = double>
auto mlmcFunction(int l, int N, const Problem& problem, const Sampler& sampler)
{
    std::array<Scalar, 6> sums = {};
    Scalar runtimeCost = 0.0;

    // TODO this loop can be parallelized
    for (int i = 0; i < N; ++i)
    {
        // sample parameters with sampler
        const auto [sFine, sCoarse] = sampler(1, l);

        // measure runtime for cost estimation
        const auto start{std::chrono::steady_clock::now()};

        // evaluate solution on both levels
        // TODO what if we have multiple quantities of interest, i.e. p is a vector?
        const auto pFine = problem.evaluate(l, sFine);
        const auto pCoarse = [&, sCoarse=sCoarse]{
            if (l == 0)
                return 0.0;
            else
                return problem.evaluate(l-1, sCoarse);
        }();

        // add runtime to total runtime cost
        const auto end{std::chrono::steady_clock::now()};
        const std::chrono::duration<Scalar> elapsedSeconds{end - start};
        runtimeCost += elapsedSeconds.count();

        // take care when parallelizing this reduction
        const auto diff = pFine-pCoarse;
        sums[0] += diff;
        sums[1] += diff*diff;
        sums[2] += diff*diff*diff;
        sums[3] += diff*diff*diff*diff;
        sums[4] += pFine;
        sums[5] += pFine*pFine;
    }

    return std::make_tuple(sums, runtimeCost);
}

/*!
 * \file
 * \ingroup UncertaintyQuantification
 * \brief MLMC estimation
 *
    Usage: (P, Nl, Cl) = mlmc(...)

    Inputs:
      N0:   initial number of samples    >  0
      eps:  desired accuracy (rms error) >  0
      Lmin: minimum level of refinement  >= 2
      Lmax: maximum level of refinement  >= Lmin

      mlmcFunction: the user low-level routine for level l estimator. Its interface is

        (sums, cost) = mlmcFunction(l, N)

        Inputs:  l: level
                 N: number of paths

        Outputs: sums[0]: sum(Y)
                 sums[1]: sum(Y**2)
                    where Y are iid samples with expected value
                        E[P_0]            on level 0
                        E[P_l - P_{l-1}]  on level l > 0
                 cost: cost of N samples

      alpha ->  weak error is  O(2^{-alpha*l})
      beta  ->  variance is    O(2^{-beta*l})
      gamma ->  sample cost is O(2^{ gamma*l})

      If alpha, beta are not positive then they will be estimated.

    Outputs:
      P:  value
      Nl: number of samples at each level
      Cl: cost of samples at each level
 */
template<class Scalar, class Functional>
auto mlmc(int Lmin, int Lmax, int N0, Scalar eps, const Functional& mlmcFunction,
          Scalar alpha0 = -1.0, Scalar beta0 = -1.0, Scalar gamma0 = -1.0)
{
    // some argument checks
    if (Lmin < 2)
        DUNE_THROW(Dumux::ParameterException, "Need Lmin >= 2");
    if (Lmax < Lmin)
        DUNE_THROW(Dumux::ParameterException, "Need Lmax >= Lmin");
    if (N0 <= 0 or eps <= 0.0)
        DUNE_THROW(Dumux::ParameterException, "Need N0 > 0, eps > 0");

    // initialization
    using std::max;
    const auto alpha = max<Scalar>(0.0, alpha0);
    const auto beta  = max<Scalar>(0.0, beta0);
    const auto gamma = max<Scalar>(0.0, gamma0);

    Scalar theta = 0.25; // what is this?
    Scalar L = Lmin;

    // allocate memory
    std::vector<int> Nl(L+1, 0.0);
    std::vector<int> dNl(L+1, N0);

    auto suml = std::array{
        std::vector<Scalar>(Lmax+1, 0.0),
        std::vector<Scalar>(Lmax+1, 0.0)
    };
    std::vector<Scalar> costl(Lmax+1, 0.0);
    std::vector<Scalar> meanl(Lmax+1, 0.0);
    std::vector<Scalar> variancel(Lmax+1, 0.0);
    std::vector<Scalar> Cl(L+1, 0.0);

    // looping while the algorithm says that
    // we still need additional samples
    while (std::accumulate(dNl.begin(), dNl.end(), 0) > 0)
    {
        // update sample sums
        for (int l = 0; l < L+1; ++l)
        {
            if (dNl[l] > 0)
            {
                const auto [sums, cost] = mlmcFunction(l, dNl[l]);
                Nl[l] += dNl[l];
                suml[0][l] += sums[0];
                suml[1][l] += sums[1];
                costl[l] += cost;
            }
        }

        // compute absolute average, variance and cost
        using std::abs, std::max;
        for (int l = 0; l < L+1; ++l)
        {
            meanl[l] = abs(suml[0][l]/Nl[l]);
            variancel[l] = max(0.0, suml[1][l]/Nl[l] - meanl[l]*meanl[l]);
            Cl[l] = costl[l]/Nl[l];
        }

        // fix to cope with possible zero values for meanl and variancel
        // (can happen in some applications when there are few samples)
        for (int l = 3; l < L+2; ++l)
        {
            meanl[l-1] = max(meanl[l-1], 0.5*meanl[l-2]/pow(2, alpha));
            variancel[l-1] = max(variancel[l-1], 0.5*variancel[l-2]/pow(2, beta));
        }

        // use linear regression to estimate alpha, beta, gamma if not given
        if (alpha0 <= 0)
            DUNE_THROW(Dune::NotImplemented, "Least squares estimater for alpha");
        if (beta0 <= 0)
            DUNE_THROW(Dune::NotImplemented, "Least squares estimater for beta");
        if (gamma0 <= 0)
            DUNE_THROW(Dune::NotImplemented, "Least squares estimater for gamma");

        // update the optimal number of additional samples dNl
        const auto updateOptimalAdditionalSamples = [&]
        {
            Scalar sumSqrtVlCl = 0.0;
            for (int l = 0; l < L+1; ++l)
                sumSqrtVlCl += sqrt(variancel[l]*Cl[l]);

            using std::ceil;
            for (int l = 0; l < L+1; ++l)
            {
                const int Ns = static_cast<int>(ceil(
                    sqrt(variancel[l]/Cl[l]) * sumSqrtVlCl / ((1.0-theta)*eps*eps)
                ));
                dNl[l] = max(0, Ns-Nl[l]);
            }
        };

        updateOptimalAdditionalSamples();

        // if (almost) converged, estimate remaining error and decide
        // whether a new level is required
        bool almostConverged = true;
        for (int l = 0; l < L+1; ++l)
        {
            if (Scalar(dNl[l]) > 0.01*Nl[l])
            {
                almostConverged = false;
                break;
            }
        }

        if (almostConverged)
        {
            // this change copes with cases with erratic meanl values
            // rang = list(range(min(3, L)))
            // rem = ( numpy.amax(ml[[L-x for x in rang]] / 2.0**(numpy.array(rang)*alpha))
            //        / (2.0**alpha - 1.0) )

            using std::pow, std::sqrt;
            const auto rem = meanl[L] / (pow(2.0, alpha) - 1.0);
            if (rem > sqrt(theta)*eps)
            {
                if (L == Lmax)
                    DUNE_THROW(Dumux::NumericalProblem, "Failed to achieve weak convergence");

                // add level and initialize values
                L++;
                meanl[L] = meanl[L-1] / pow(2.0, alpha);
                variancel[L] = variancel[L-1] / pow(2.0, beta);

                Cl.push_back(Cl.back() * pow(2.0, gamma));
                Nl.push_back(0);
                dNl.push_back(0);

                updateOptimalAdditionalSamples();
            }
        }
    }

    // finally, evaluate the multi-level estimator
    // TODO what if more than 1 quantitiy of interest?
    // TODO how to add a multi-level estimator for the variance?
    Scalar P = 0.0;
    for (int l = 0; l < L+1; ++l)
        P += suml[0][l]/Nl[l];

    return std::make_tuple(P, Nl, Cl);
}

} // end namespace Dumux

#endif
