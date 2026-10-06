// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Fučík's semi-analytical McWhorter-Sunada solution for counter-current imbibition.
 */
#ifndef DUMUX_TEST_TWOP_MCWHORTERSUNADA_ANALYTICSOLUTION_HH
#define DUMUX_TEST_TWOP_MCWHORTERSUNADA_ANALYTICSOLUTION_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>

namespace Dumux {

/*!
 * \brief Saturation reference for the homogeneous, gravity-free, incompressible test.
 *
 * Implements method B of Fučík et al. @cite Fucik2007, specialized to R = 0
 * (closed right boundary). Here, the notation used in @cite Fucik2007 is adopted.
 * The default benchmark end time is chosen before the semi-infinite reference front reaches the right boundary.
 * See README.md for equations, quadrature and benchmark assumptions.
 */
template<class TypeTag>
class McWhorterAnalyticSolution
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

public:
    explicit McWhorterAnalyticSolution(std::shared_ptr<const Problem> problem,
                                      int intervals = getParam<int>("Reference.NumIntervals", 4000))
    : problem_(std::move(problem))
    , values_(problem_->gridGeometry().numDofs())
    {
        static_assert(ModelTraits::priVarFormulation() == TwoPFormulation::p0s1,
                      "This benchmark uses wetting pressure and nonwetting saturation");
        if (getParam<bool>("Problem.EnableGravity", true))
            DUNE_THROW(Dune::InvalidStateException, "McWhorter reference requires gravity to be disabled");

        const auto pos = problem_->gridGeometry().bBoxMin();
        const auto& spatialParams = problem_->spatialParams();
        const auto interaction = spatialParams.fluidMatrixInteractionAtPos(pos);
        const auto& residual = interaction.pcSwCurve().effToAbsParams();
        swInitial_ = residual.swr(); // paper notation: S_i
        swBoundary_ = 1.0 - residual.snr(); // paper notation: S_0
        xMin_ = pos[0];
        porosity_ = spatialParams.porosityAtPos(pos);
        const Scalar permeability = spatialParams.permeabilityAtPos(pos);
        const Scalar tolerance = getParam<Scalar>("Reference.Tolerance", 1e-10);
        const int maxIterations = getParam<int>("Reference.MaxIterations", 10000);
        if (intervals < 2 || !(tolerance > 0.0) || maxIterations < 1 || !(swBoundary_ > swInitial_) || !(porosity_ > 0.0) || !(permeability > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "Invalid McWhorter reference parameters");

        FluidState fluidState;
        fluidState.setTemperature(spatialParams.temperatureAtPos(pos));
        fluidState.setPressure(FluidSystem::phase0Idx, problem_->referencePressure());
        fluidState.setPressure(FluidSystem::phase1Idx, problem_->referencePressure());
        const Scalar muW = FluidSystem::viscosity(fluidState, FluidSystem::phase0Idx);
        const Scalar muN = FluidSystem::viscosity(fluidState, FluidSystem::phase1Idx);
        if (!(muW > 0.0) || !(muN > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "McWhorter reference requires positive viscosities");

        const Scalar h = (swBoundary_ - swInitial_)/intervals;
        saturation_.resize(intervals + 1);
        xi_.resize(intervals + 1);
        std::vector<Scalar> D(intervals + 1); // D(S): capillary diffusivity
        std::vector<Scalar> G(intervals + 1); // G(S): integral-equation integrand
        std::vector<Scalar> cumulativeWeightedIntegral(intervals + 1);
        std::vector<Scalar> cumulativeTailIntegral(intervals + 1);
        // Fill D(S)
        for (int i = 0; i <= intervals; ++i)
        {
            const Scalar sw = swInitial_ + i*h;
            saturation_[i] = sw;
            const Scalar lambdaW = interaction.krw(sw)/muW;
            const Scalar lambdaN = interaction.krn(sw)/muN;
            D[i] = -permeability*lambdaW*lambdaN/(lambdaW + lambdaN)*interaction.dpc_dsw(sw);
            if (!std::isfinite(D[i]) || D[i] < 0.0)
                DUNE_THROW(Dune::InvalidStateException, "Invalid capillary diffusivity in McWhorter reference");
        }
        // Both mobilities vanish at their respective residual endpoints. The limit
        // G(S_i) is zero for the regularized Brooks-Corey law used by this test.
        D.front() = 0.0;
        D.back() = 0.0;
        G = D; // G_0 = D because F_0 = 1 for R = 0

        // Integrate the piecewise-linear G(S). The weighted integral is
        // integral_{S_i}^{S} (s - S_i) G(s) ds and the tail integral is
        // integral_{S}^{S_0} G(s) ds.
        auto computeCumulativeIntegrals = [&]
        {
            cumulativeWeightedIntegral[0] = 0.0;
            cumulativeTailIntegral[intervals] = 0.0;
            for (int i = 0; i < intervals; ++i)
                cumulativeWeightedIntegral[i+1] = cumulativeWeightedIntegral[i] + h*h/6.0*((3*i + 1)*G[i] + (3*i + 2)*G[i+1]);
            for (int i = intervals - 1; i >= 0; --i)
                cumulativeTailIntegral[i] = cumulativeTailIntegral[i+1] + 0.5*h*(G[i] + G[i+1]);
        };

        bool converged = false;
        for (int iteration = 0; iteration < maxIterations; ++iteration)
        {
            computeCumulativeIntegrals();
            const Scalar integral = cumulativeWeightedIntegral.back(); // I(S_0), normalization integral
            if (!(integral > 0.0) || !std::isfinite(integral))
                DUNE_THROW(Dune::InvalidStateException, "Degenerate McWhorter reference integral");
            Scalar change = 0.0, scale = 0.0;
            for (int i = 1; i < intervals; ++i)
            {
                // F(S): normalized profile function, evaluated without cancellation near S_i.
                const Scalar F = (cumulativeWeightedIntegral[i] + i*h*cumulativeTailIntegral[i])/integral;
                const Scalar next = D[i]/F; // G = D/F for Fučík's R = 0 specialization
                if (!(F > 0.0) || !std::isfinite(next))
                    DUNE_THROW(Dune::InvalidStateException, "Invalid Fučík iteration");
                change = std::max(change, std::abs(next - G[i]));
                scale = std::max(scale, std::abs(next));
                G[i] = next;
            }
            if (change <= tolerance*scale)
            {
                converged = true;
                break;
            }
        }
        if (!converged)
            DUNE_THROW(Dune::InvalidStateException, "Fučík reference iteration did not converge");

        computeCumulativeIntegrals();
        const Scalar integral = cumulativeWeightedIntegral.back(); // I(S_0), normalization integral
        A_ = std::sqrt(0.5*porosity_*integral); // A: amplitude in the inverted profile
        for (int i = 0; i <= intervals; ++i)
            xi_[i] = 2.0*A_/porosity_*cumulativeTailIntegral[i]/integral;
        update(0.0);
    }

    Scalar computeSaturation(Scalar x, Scalar time) const
    {
        if (time <= 0.0)
            return swInitial_;
        const Scalar xi = (x - xMin_)/std::sqrt(time);
        if (xi <= 0.0)
            return swBoundary_;
        if (xi >= xi_.front())
            return swInitial_;
        const auto upper = std::lower_bound(xi_.begin(), xi_.end(), xi, std::greater<Scalar>());
        const auto i = std::distance(xi_.begin(), upper);
        const Scalar weight = (xi_[i-1] - xi)/(xi_[i-1] - xi_[i]);
        return saturation_[i-1] + weight*(saturation_[i] - saturation_[i-1]);
    }

    void update(Scalar time)
    {
        for (const auto& element : elements(problem_->gridGeometry().gridView()))
            values_[problem_->gridGeometry().elementMapper().index(element)] = computeSaturation(element.geometry().center()[0], time);
    }

    const std::vector<Scalar>& values() const { return values_; }

    Scalar initialSaturation() const { return swInitial_; }

private:
    std::shared_ptr<const Problem> problem_;
    std::vector<Scalar> values_, saturation_, xi_;
    Scalar swInitial_, swBoundary_, porosity_, xMin_, A_;
};

} // namespace Dumux
#endif
