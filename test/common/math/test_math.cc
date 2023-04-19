// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief This file tests several math functions:
 *     the vtmv function with vtmv(Vector, FieldScalar, Vector) and vtmv(Vector, Matrix, Vector).
 *     the trace function
 *     the harmonicMean function
 * \todo test more math functions!
 *
 * We declare some vectors and matrices and test the combinations of them
 * against a previously calculated result.
 */
#include <config.h>

#include <iostream>
#include <utility>
#include <algorithm>
#include <array>
#include <random>
#include <tuple>

#include <dune/common/float_cmp.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

#include <dumux/common/math.hh>
#include <dumux/io/format.hh>

namespace Test {

template<class Scalar, class Table>
void checkTableInterpolation(Scalar ip, Scalar expected, const Table& table, Scalar eps = 1e-15)
{
    using namespace Dumux;
    auto interp = interpolate<InterpolationPolicy::LinearTable>(ip, table);
    if (!Dune::FloatCmp::eq(interp, expected, eps))
        DUNE_THROW(Dune::Exception, "Wrong interpolation, expected " << expected << " got " << interp);
}

} // end namespace Test

int main()
{
    //! Declare Vectors (FieldVector and DynamicVector)
    Dune::FieldVector<double, 3> v1({1.0, 2.0, 3.0});
    Dune::FieldVector<double, 3> v2({1.0, 0.0, -2.0});

    Dune::DynamicVector<double> v1_dyn(v1);
    Dune::DynamicVector<double> v2_dyn(v2);

    //! Declare 3x3 Matrices with 3's on the principal diagonal
    double k = 3.0;
    Dune::FieldMatrix<double, 3, 3> K(0.0);
    K[0][0] = k;
    K[1][1] = k;
    K[2][2] = k;
    Dune::DynamicMatrix<double> K_dyn(K);


    //////////////////////////////////////////////////////////////////
    ///// Dumux::vtmv ////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////

    //! Set reference result. Should be -15 for all combinations
    const double reference = -15;

    //! Test with FieldVector v1
    //! Test vtmv function with FieldVector v1 and Scalar k and FieldVector v2
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1, k, v2), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");
    //! Test vtmv function with FieldVector v1 and FieldMatrix K and FieldVector v2
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1, K, v2), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");
    //! Test vtmv function with FieldVector v1 and FieldMatrix K and DynamicVector v2_dyn
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1, K, v2_dyn), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");
    //! Test vtmv function with FieldVector v1 and DynamicMatrix K_dyn and FieldVector v2
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1, K_dyn, v2), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");

    //! Test with DynamicVector v1_dyn
    //! Test vtmv function with DynamicVector v1_dyn and Scalar k and FieldVector v2
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1_dyn, k, v2), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");
    //! Test vtmv function with DynamicVector v1_dyn and FieldMatrix K and FieldVector v2
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1_dyn, K, v2), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");
    //! Test vtmv function with DynamicVector v1_dyn and Scalar k and DynamicVector v2_dyn
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1_dyn, K, v2_dyn), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");
    //! Test vtmv function with DynamicVector v1_dyn and DynamicMatrix K_dyn and FieldVector v2
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1_dyn, K_dyn, v2), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");

    //! Test Dynamic Vectors and Dynamic Matrices
    //! Test vtmv function with DynamicVector v1_dyn and DynamicMatrix K_dyn and DynamicVector v2_dyn
    if (!Dune::FloatCmp::eq(reference, Dumux::vtmv(v1_dyn, K_dyn, v2_dyn), 1e-6))
        DUNE_THROW(Dune::Exception, "vtmv-result does not match reference");

    //////////////////////////////////////////////////////////////////
    ///// Dumux::mv ///////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////

    const Dune::FieldVector<double, 3> fieldVecRef({3.0, 6.0, 9.0});
    const Dune::DynamicVector<double> dynVecRef({3.0, 6.0, 9.0});

    //! Test mv with FieldMatrix and FieldVector
    if (!Dune::FloatCmp::eq<Dune::FieldVector<double, 3>, Dune::FloatCmp::CmpStyle::absolute>(fieldVecRef, Dumux::mv(K, v1), 1e-6))
        DUNE_THROW(Dune::Exception, "mv-result does not match reference");

    //! Test mv with Scalar and FieldVector
    if (!Dune::FloatCmp::eq<Dune::FieldVector<double, 3>, Dune::FloatCmp::CmpStyle::absolute>(fieldVecRef, Dumux::mv(k, v1), 1e-6))
        DUNE_THROW(Dune::Exception, "mv-result does not match reference");

    //! Test mv with FieldMatrix and DynamicVector (use cast to FieldVector, there is no specialization of eq() for DynamicVector)
    if (!Dune::FloatCmp::eq<Dune::FieldVector<double, 3>, Dune::FloatCmp::CmpStyle::absolute>(dynVecRef, Dumux::mv(K, v1_dyn), 1e-6))
        DUNE_THROW(Dune::Exception, "mv-result does not match reference");

    //! Test mv with Scalar and DynamicVector (use cast to FieldVector, there is no specialization of eq() for DynamicVector)
    if (!Dune::FloatCmp::eq<Dune::FieldVector<double, 3>, Dune::FloatCmp::CmpStyle::absolute>(dynVecRef, Dumux::mv(k, v1_dyn), 1e-6))
        DUNE_THROW(Dune::Exception, "mv-result does not match reference");

    //! Test mv with DynamicMatrix and DynamicVector (use cast to FieldVector, there is no specialization of eq() for DynamicVector)
    if (!Dune::FloatCmp::eq<Dune::FieldVector<double, 3>, Dune::FloatCmp::CmpStyle::absolute>(dynVecRef, Dumux::mv(K_dyn, v1_dyn), 1e-6))
        DUNE_THROW(Dune::Exception, "mv-result does not match reference");

    //! Test mv with DynamicMatrix and FieldVector
    if (!Dune::FloatCmp::eq<Dune::FieldVector<double, 3>, Dune::FloatCmp::CmpStyle::absolute>(fieldVecRef, Dumux::mv(K_dyn, v1), 1e-6))
        DUNE_THROW(Dune::Exception, "mv-result does not match reference");

    //////////////////////////////////////////////////////////////////
    ///// Dumux::trace ///////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////

    const auto t1 = Dumux::trace(K);
    const auto t2 = Dumux::trace(K_dyn);
    if (!Dune::FloatCmp::eq(t1, t2, 1e-30))
        DUNE_THROW(Dune::Exception, "Traces do not match!");


    //////////////////////////////////////////////////////////////////
    ///// Dumux::harmonicMean ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////

    constexpr auto mean = Dumux::harmonicMean(4.0, 5.0);
    static_assert( (mean - 4.0*5.0*2.0/(4.0+5.0)) < 1e-30 && (mean - 4.0*5.0*2.0/(4.0+5.0)) > -1e-30 , "Wrong harmonic mean!");

    //////////////////////////////////////////////////////////////////
    ///// Dumux::sign ////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////

    static_assert(Dumux::sign(0.0) == 0, "Wrong sign!");
    static_assert(Dumux::sign(-0.0) == 0, "Wrong sign!");
    static_assert(Dumux::sign(0) == 0, "Wrong sign!");
    static_assert(Dumux::sign(-0) == 0, "Wrong sign!");
    static_assert(Dumux::sign(1) == 1, "Wrong sign!");
    static_assert(Dumux::sign(2.0) == 1, "Wrong sign!");
    static_assert(Dumux::sign(-2) == -1, "Wrong sign!");
    static_assert(Dumux::sign(-3.5) == -1, "Wrong sign!");

    //////////////////////////////////////////////////////////////////
    ///// Dumux::interpolate /////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    std::vector<double> a{0.0, 1.0, 2.0};
    std::vector<double> b{-1.0, 1.0, 3.0};
    const auto table = std::make_pair(a, b);

    Test::checkTableInterpolation(-1.0, -1.0, table);
    Test::checkTableInterpolation(+0.0, -1.0, table);
    Test::checkTableInterpolation(0.001, -0.998, table);
    Test::checkTableInterpolation(1.5, 2.0, table);
    Test::checkTableInterpolation(2.0, 3.0, table);
    Test::checkTableInterpolation(3.0, 3.0, table);

    //////////////////////////////////////////////////////////////////
    ///// Dumux::linspace ////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    {
        const auto v = Dumux::linspace(1.0, 2.0, 100);
        if (!Dune::FloatCmp::eq(v.back(), 2.0))
            DUNE_THROW(Dune::Exception, "[linspace] Last point not included!");
        if (!Dune::FloatCmp::eq(v.front(), 1.0))
            DUNE_THROW(Dune::Exception, "[linspace] First point not included!");
        if (v.size() != 100)
            DUNE_THROW(Dune::Exception, "[linspace] Size incorrect!");
        if (!std::is_sorted(v.begin(), v.end()))
            DUNE_THROW(Dune::Exception, "[linspace] Not sorted in ascending order!");
    }
    {
        const auto v = Dumux::linspace(1.0, 2.0, 100, /*endPoint=*/false);
        if (!Dune::FloatCmp::eq(v.back(), 2.0-0.01))
            DUNE_THROW(Dune::Exception, "[linspace] Last point not correct!");
        if (!Dune::FloatCmp::eq(v.front(), 1.0))
            DUNE_THROW(Dune::Exception, "[linspace] First point not included!");
        if (v.size() != 100)
            DUNE_THROW(Dune::Exception, "[linspace] Size incorrect!");
        if (!std::is_sorted(v.begin(), v.end()))
            DUNE_THROW(Dune::Exception, "[linspace] Not sorted in ascending order!");
    }

    //////////////////////////////////////////////////////////////////
    ///// Dumux::linearRegression ////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    {
        const auto x = Dumux::linspace(0.0, 1.0, 100);
        const auto [intercept, slope] = Dumux::linearRegression(x, x);
        if (std::abs(intercept) > 1e-14)
            DUNE_THROW(Dune::Exception, "[linearRegreesion] Wrong intercept " << intercept << ", should be 0.0.");
        if (std::abs(slope-1.0) > 1e-14)
            DUNE_THROW(Dune::Exception, "[linearRegreesion] Wrong slope " << slope << ", should be 1.0.");
    }
    {
        const auto x = Dumux::linspace(0.0, 1.0, 100);
        const auto y = Dumux::linspace(1.0, 3.0, 100);
        const auto [intercept, slope] = Dumux::linearRegression(x, y);
        if (std::abs(intercept-1.0) > 1e-14)
            DUNE_THROW(Dune::Exception, "[linearRegreesion] Wrong intercept " << intercept << ", should be 1.0.");
        if (std::abs(slope-2.0) > 1e-14)
            DUNE_THROW(Dune::Exception, "[linearRegreesion] Wrong slope " << slope << ", should be 2.0.");
    }

    //////////////////////////////////////////////////////////////////
    ///// Dumux::invertCubicPolynomial ///////////////////////////////
    //////////////////////////////////////////////////////////////////

    const auto cubicCoefficientsFromRoots = [](auto x, auto y, auto z, double scale = 1.0){
        const auto b = -(x + y + z);
        const auto c = x*y + x*z + y*z;
        const auto d = -x*y*z;
        // scaling doesn't change the roots
        return std::make_tuple(1.0*scale, b*scale, c*scale, d*scale);
    };

    const auto quadCoefficientsFromRoots = [](auto x, auto y, double scale = 1.0){
        const auto b = -(x + y);
        const auto c = x*y;
        return std::make_tuple(1.0*scale, b*scale, c*scale);
    };

    const auto linearCoefficientsFromRoots = [](auto x, double scale = 1.0){
        return std::make_tuple(1.0*scale, -x*scale);
    };

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(-10.0, 10.0);
    auto threeRandomNumbers = [&]() mutable {
        const auto x = dis(gen);
        auto y = dis(gen); while (Dune::FloatCmp::eq(x, y, 1e-3)) { y = dis(gen); }
        auto z = dis(gen); while (Dune::FloatCmp::eq(x, z, 1e-3) || Dune::FloatCmp::eq(y, z, 1e-3)) { z = dis(gen); }
        return std::make_tuple(x, y, z);
    };

    const auto checkRoots = [](const std::array<double, 3>& roots,
                               const std::array<double, 3>& ref,
                               const std::array<double, 4>& coeff,
                               int n, bool throwing = true) -> bool
    {
        const auto refToString = [&](){
            std::string refStr{};
            for (int i = 0; i < n; i++)
                refStr += Dumux::Fmt::format("{:.14e},", ref[i]);
            return refStr;
        };

        for (int i = 0; i < n; i++)
        {
            using std::isfinite;
            if (!isfinite(roots[i]))
                DUNE_THROW(Dune::Exception, "Expecting finite root");
            if (throwing && std::none_of(ref.begin(), ref.end(), [&](auto r){ return Dune::FloatCmp::eq(roots[i], r, 1e-5); }))
                DUNE_THROW(Dune::Exception, "[invertCubicPolynomial] Root " << Dumux::Fmt::format(
                        "{} of {}: {:.14e}, does not match reference roots = [{}].\n", i+1, n, roots[i], refToString())
                        << Dumux::Fmt::format("Polynomial: {:+.5e}x³ {:+.5e}x² {:+.5e}x {:+.5e}\n", coeff[0], coeff[1], coeff[2], coeff[3]));
            else
                return false;
        }

        return true;
    };

    std::array<double, 3> roots{};
    for (double scaling : {1.0, -1.0, 5.0})
    {
        // random numbers
        for (int i = 0; i < 100000; ++i)
        {
            const auto [x, y, z] = threeRandomNumbers();

            // test for three unique roots
            {
                const auto [a, b, c, d] = cubicCoefficientsFromRoots(x, y, z, scaling);
                int numRoots = Dumux::invertCubicPolynomial(roots.data(), a, b, c, d);
                if (numRoots != 3)
                    DUNE_THROW(Dune::Exception, "Expecting three roots! Got " << numRoots);
                if (!checkRoots(roots, {x, y, z}, {a, b, c, d}, numRoots, false))
                {
                    // Try to increase Newton iterations for increased precision and try again
                    int numRoots = Dumux::invertCubicPolynomial(roots.data(), a, b, c, d, 10);
                    checkRoots(roots, {x, y, z}, {a, b, c, d}, numRoots);
                }
            }

            // test for two unique roots
            {
                const auto [b, c, d] = quadCoefficientsFromRoots(x, z, scaling);
                int numRoots = Dumux::invertCubicPolynomial(roots.data(), 0.0, b, c, d);
                if (numRoots != 2)
                    DUNE_THROW(Dune::Exception, "Expecting two roots! Got " << numRoots);
                if (!checkRoots(roots, {x, z}, {0.0, b, c, d}, numRoots, false))
                {
                    // Try to increase Newton iterations for increased precision and try again
                    int numRoots = Dumux::invertCubicPolynomial(roots.data(), 0.0, b, c, d, 10);
                    checkRoots(roots, {x, z}, {0.0, b, c, d}, numRoots);
                }
            }

            // test for one unique root
            {
                const auto [c, d] = linearCoefficientsFromRoots(x, scaling);
                int numRoots = Dumux::invertCubicPolynomial(roots.data(), 0.0, 0.0, c, d);
                if (numRoots != 1)
                    DUNE_THROW(Dune::Exception, "Expecting one root! Got " << numRoots);
                if (!checkRoots(roots, {x}, {0.0, 0.0, c, d}, numRoots, false))
                {
                    // Try to increase Newton iterations for increased precision and try again
                    int numRoots = Dumux::invertCubicPolynomial(roots.data(), 0.0, 0.0, c, d, 10);
                    checkRoots(roots, {x}, {0.0, 0.0, c, d}, numRoots);
                }
            }
        }

        // some corner cases with close roots
        {
            const auto [x, y, z] = std::make_tuple(10.0, 10.0+1e-8, 0.0);
            const auto [a, b, c, d] = cubicCoefficientsFromRoots(x, y, z, scaling);
            int numRoots = Dumux::invertCubicPolynomial(roots.data(), a, b, c, d);
            if (numRoots != 3)
                DUNE_THROW(Dune::Exception, "Expecting three roots! Got " << numRoots);
            if (!checkRoots(roots, {x, y, z}, {a, b, c, d}, numRoots, false))
            {
                // Try to increase Newton iterations for increased precision and try again
                int numRoots = Dumux::invertCubicPolynomial(roots.data(), a, b, c, d, 10);
                checkRoots(roots, {x, y, z}, {a, b, c, d}, numRoots);
            }
        }{
            const auto [x, y, z] = std::make_tuple(10.0, 10.0+1e-4, 10.0+2e-4);
            const auto [a, b, c, d] = cubicCoefficientsFromRoots(x, y, z, scaling);
            int numRoots = Dumux::invertCubicPolynomial(roots.data(), a, b, c, d);
            if (numRoots != 3)
                DUNE_THROW(Dune::Exception, "Expecting three roots! Got " << numRoots);
            if (!checkRoots(roots, {x, y, z}, {a, b, c, d}, numRoots, false))
            {
                // Try to increase Newton iterations for increased precision and try again
                int numRoots = Dumux::invertCubicPolynomial(roots.data(), a, b, c, d, 10);
                checkRoots(roots, {x, y, z}, {a, b, c, d}, numRoots);
            }
        }
    }

    // Test one corner case (see MR!2594)
    int nroots = Dumux::invertCubicPolynomial(roots.data(),
        1.0, -0.96165703943410097, 0.30826068470787077, 0.01340255221155587
    );
    if (nroots != 1)
        DUNE_THROW(Dune::Exception, "Expecting one root");
    using std::isfinite;
    if (!isfinite(roots[0]))
        DUNE_THROW(Dune::Exception, "Expecting finite root");
    if (!Dune::FloatCmp::eq(roots[0], -3.863448718244389e-02, 1e-14))
        DUNE_THROW(Dune::Exception, "Root of cubic equation does not match reference");

    return 0;
}
