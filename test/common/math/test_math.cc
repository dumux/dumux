// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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

#include <dune/common/float_cmp.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

#include <dumux/common/math.hh>

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

}
