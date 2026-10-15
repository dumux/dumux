// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Linear stability of the backward differentiation formulas (BDF)
 *
 * For the Dahlquist test equation du/dt = lambda u, i.e. M(u) = u and R(u) = -lambda u, a
 * k-step method with constant step size gives the recurrence sum_j (alpha_j - z beta_j) u^{n+1-j} = 0
 * with z = lambda dt. The method is stable for z if all roots of
 * p(zeta) = sum_j (alpha_j - z beta_j) zeta^{k-j} lie in the closed unit disk. The boundary of
 * the stability region is the root locus z(theta) = sum_j alpha_j exp(-i j theta) (beta = (1, 0, ...)).
 *
 * The test checks that the roots on the locus have modulus one, A-stability of BDF1 and BDF2,
 * the A(alpha) angles of BDF3 to BDF5 against the values in Hairer & Wanner (1996), the damping
 * of stiff modes, and the parasitic root omega^2/(1 + 2 omega) of BDF2 for a constant step size
 * ratio omega. It writes the data for the stability figure.
 */
#include <config.h>

#include <cmath>
#include <array>
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include <dune/common/exceptions.hh>

#include <dumux/io/format.hh>
#include <dumux/io/json.hh>
#include <dumux/common/initialize.hh>
#include <dumux/experimental/timestepping/multistepmethods.hh>

#include "polynomialroots.hh"

namespace Dumux {

using Scalar = double;
using Complex = std::complex<Scalar>;

void expect(const bool condition, const std::string& message)
{
    if (!condition)
        DUNE_THROW(Dune::InvalidStateException, message);
}

void expectNear(const Scalar value, const Scalar reference, const Scalar tolerance, const std::string& message)
{
    using std::abs;
    if (!(abs(value - reference) <= tolerance))
        DUNE_THROW(Dune::InvalidStateException, message << " (" << value << " vs. " << reference << ")");
}

//! The largest modulus of the roots of the characteristic polynomial for z = lambda dt
Scalar rootModulus(const std::vector<Scalar>& alpha, const Complex z)
{
    std::vector<Complex> coefficients(alpha.begin(), alpha.end());
    coefficients[0] -= z;
    Scalar modulus = 0.0;
    for (const auto& r : Test::roots(coefficients))
        modulus = std::max(modulus, std::abs(r));
    return modulus;
}

//! The boundary locus z(theta) = sum_j alpha_j exp(-i j theta)
std::vector<Complex> boundaryLocus(const std::vector<Scalar>& alpha, const int numPoints = 2001)
{
    std::vector<Complex> locus(numPoints);
    for (int i = 0; i < numPoints; ++i)
    {
        const auto theta = 2.0*M_PI*i/(numPoints - 1);
        Complex z = 0.0;
        for (std::size_t j = 0; j < alpha.size(); ++j)
            z += alpha[j]*std::polar(1.0, -Scalar(j)*theta);
        locus[i] = z;
    }
    return locus;
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    Dumux::initialize(argc, argv);

    using BDF = Experimental::MultiStep::BDF<Scalar>;

    // A(alpha) angles in degrees, Hairer & Wanner (1996), Section V.2
    const std::array<Scalar, 5> referenceAngles = {90.0, 90.0, 86.03, 73.35, 51.84};

    Dumux::Json::JsonTree out;
    out["methods"] = Dumux::Json::JsonTree::object();

    for (std::size_t k = 1; k <= 5; ++k)
    {
        const BDF bdf(k);
        std::vector<Scalar> alpha(k+1), beta(k+1);
        bdf.weights(std::vector<Scalar>(k, 1.0), alpha, beta);

        const auto locus = boundaryLocus(alpha);
        Scalar maxDeviation = 0.0;
        for (std::size_t i = 0; i < locus.size(); i += 20)
            maxDeviation = std::max(maxDeviation, std::abs(rootModulus(alpha, locus[i]) - 1.0));
        expect(maxDeviation < 1e-8, Fmt::format("The boundary locus of BDF{} has roots off the unit circle", k));

        // the A(alpha) angle is the smallest angle of a locus point in the left half plane to the negative real axis
        Scalar angle = 90.0;
        for (const auto& z : locus)
            if (z.real() < -1e-12)
                angle = std::min(angle, std::abs(std::arg(-z))*180.0/M_PI);
        std::cout << Fmt::format("-- BDF{}: A(alpha) angle {:.2f} degrees\n", k, angle);
        expectNear(angle, referenceAngles[k-1], 0.01, Fmt::format("A(alpha) angle of BDF{}", k));

        // stable inside the sector |arg(-z)| < alpha (sampled up to |z| = 1e4)
        Scalar maxModulusInSector = 0.0;
        for (int i = 0; i <= 40; ++i)
        {
            const auto radius = std::pow(10.0, -2.0 + 6.0*i/40.0);
            for (int j = 0; j <= 40; ++j)
            {
                const auto phi = (angle - 0.5)*M_PI/180.0*(-1.0 + 2.0*j/40.0);
                maxModulusInSector = std::max(maxModulusInSector, rootModulus(alpha, -std::polar(radius, phi)));
            }
        }
        expect(maxModulusInSector <= 1.0 + 1e-10, Fmt::format("BDF{} is unstable inside its A(alpha) sector", k));

        // stiff modes are damped
        const auto stiffModulus = rootModulus(alpha, Complex{-1e6, 0.0});
        std::cout << Fmt::format("-- BDF{}: largest root modulus for lambda dt = -1e6: {:.3e}\n", k, stiffModulus);
        expect(stiffModulus < 0.1, Fmt::format("BDF{} should damp stiff modes", k));

        auto& node = out["methods"][Fmt::format("bdf{}", k)];
        node["name"] = bdf.name();
        node["angle"] = angle;
        std::vector<std::array<Scalar, 2>> points;
        for (const auto& z : locus)
            points.push_back({z.real(), z.imag()});
        node["boundaryLocus"] = points;
    }

    // parasitic root of BDF2 for the constant step size ratio omega (u' = 0)
    std::vector<std::array<Scalar, 2>> parasitic;
    for (int i = 0; i <= 300; ++i)
    {
        const Scalar omega = 0.05 + 2.95*i/300.0;
        const BDF bdf(2);
        std::vector<Scalar> alpha(3), beta(3);
        bdf.weights(std::vector<Scalar>{omega, 1.0}, alpha, beta);
        const auto r = Test::roots({alpha[0], alpha[1], alpha[2]});
        const auto modulus = std::abs(std::abs(r[0] - 1.0) < std::abs(r[1] - 1.0) ? r[1] : r[0]);
        expectNear(modulus, omega*omega/(1.0 + 2.0*omega), 1e-10, "Parasitic root of variable-step BDF2");
        parasitic.push_back({omega, modulus});
    }
    out["bdf2ParasiticRoot"] = parasitic;
    out["bdf2MaxStepSizeRatioDefault"] = 2.0;

    const std::string outputFile = "test_multistepmethods_stability_data.json";
    std::ofstream output(outputFile);
    output << std::setw(2) << out << std::endl;
    std::cout << "Wrote stability data to " << outputFile << std::endl;

    return 0;
}
