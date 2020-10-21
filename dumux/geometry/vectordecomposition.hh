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

#ifndef DUMUX_GEOMETRY_VECTORDECOMPOSITION_HH
#define DUMUX_GEOMETRY_VECTORDECOMPOSITION_HH

#include <vector>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dumux/common/math.hh>

namespace Dumux {

template<class V>
std::pair<bool, std::array<typename V::value_type, V::dimension>>
expressInBasis(const V& a, const V& b)
{
    std::array<typename V::value_type, V::dimension> coeff;
    coeff.fill(0.0);

    // if they are colinear this "decomposition" works
    const auto aNorm2 = a.two_norm2();
    const auto bNorm2 = b.two_norm2();
    const auto ab = a*b;
    using std::max;
    const auto eps2 = 1e-14*max(aNorm2, bNorm2);
    if (ab*ab > aNorm2*bNorm2 - eps2*eps2)
    {
        coeff[0] = ab/bNorm2;
        return {true, std::move(coeff)};
    }
    else
        return {false, std::move(coeff)};
}

template<class V, std::enable_if_t<V::dimension == 3, int> = 0>
std::pair<bool, std::array<typename V::value_type, V::dimension>>
expressInBasis(V a, V b1, V b2, V b3)
{
    std::array<typename V::value_type, V::dimension> coeff;
    coeff.fill(0.0);

    const auto b1Norm = b1.two_norm();
    const auto b2Norm = b2.two_norm();
    const auto b3Norm = b3.two_norm();
    const auto aNorm = a.two_norm();

    // normalize all vectors for stability
    b1 /= b1Norm;
    b2 /= b2Norm;
    b3 /= b3Norm;
    a /= aNorm;

    const auto A_f = b1*crossProduct(b2, b3);
    const auto A_f_1 = a*crossProduct(b2, b3);
    const auto A_f_2 = b1*crossProduct(a, b3);
    const auto A_f_3 = b1*crossProduct(b2, a);

    if (std::abs(A_f) < 1.0e-8)
        return {false, std::move(coeff)};

    coeff[0] = A_f_1 / A_f;
    coeff[1] = A_f_2 / A_f;
    coeff[2] = A_f_3 / A_f;

    std::transform(coeff.begin(), coeff.end(), coeff.begin(), [](const auto& c){
        using std::abs;
        return abs(c) < 1.0e-12 ? 0.0 : c;
    });

    coeff[0] *= aNorm/b1Norm;
    coeff[1] *= aNorm/b2Norm;
    coeff[2] *= aNorm/b3Norm;

    return {true, coeff};
}

template<class V>
std::tuple<std::array<int, V::dimension>, std::array<typename V::value_type, V::dimension>, std::size_t>
decomposition(const V& vec, const std::vector<V>& b, int includeIndex = -1)
{
    assert(b.size() > 0);

    using Scalar = typename V::value_type;
    static constexpr auto dim = V::dimension;
    std::array<int, dim> indices; indices.fill(0);

    {
        for (int i = 0; i < b.size(); ++i)
        {
            const auto [valid, coeff] = expressInBasis(vec, b[i]);
            if (valid && coeff[0] > 0.0)
            {
                indices = {i, 0, 0};
                return { std::move(indices), std::move(coeff), 1 };
            }
        }
    }

    if constexpr (dim == 3)
    {
        const auto numVectors = b.size();
        if (numVectors < 3)
            DUNE_THROW(Dune::NotImplemented, "Decomposing vector in 3D with less than 3 basis candidates");

        bool foundValidPositiveSolution = false;
        bool foundValidNeg1Solution = false;
        bool foundValidNeg2Solution = false;

        std::array<Scalar, V::dimension> bestCoefficients;

        std::array<int, dim> indicesNeg; indicesNeg.fill(0);
        std::array<Scalar, V::dimension> bestCoefficientsNeg;

        std::array<int, dim> indicesNeg2; indicesNeg2.fill(0);
        std::array<Scalar, V::dimension> bestCoefficientsNeg2;

        auto cost = std::numeric_limits<Scalar>::max();
        auto costNeg = std::numeric_limits<Scalar>::max();
        auto costNeg2 = std::numeric_limits<Scalar>::max();

        for (int i = 0; i < numVectors - 2; ++i)
        {
            for (int j = i+1; j < numVectors - 1; ++j)
            {
                for (int k = j+1; k < numVectors; ++k)
                {
                    // exclude combinations that don't contain the includeIndex vector
                    if (includeIndex >= 0 && !(i == includeIndex || j == includeIndex || k == includeIndex))
                        continue;

                    const auto [valid, coeff] = expressInBasis(vec, b[i], b[j], b[k]);
                    if (valid && coeff[includeIndex] > 0.0)
                    {
                        using std::signbit;
                        const auto numNeg = std::count_if(coeff.begin(), coeff.end(), [](const auto& c){ return signbit(c); });
                        if (numNeg == 0)
                        {
                            const auto coeffSum = std::accumulate(coeff.begin(), coeff.end(), 0.0);
                            foundValidPositiveSolution = true;
                            if (coeffSum < cost)
                            {
                                cost = coeffSum;
                                bestCoefficients = std::move(coeff);
                                indices = {i, j, k};
                            }
                        }
                        else if (numNeg == 1)
                        {
                            using std::abs;
                            const auto thisCost = abs(*std::min_element(coeff.begin(), coeff.end()));
                            foundValidNeg1Solution = true;

                            if (thisCost < costNeg)
                            {
                                costNeg = thisCost;
                                bestCoefficientsNeg = std::move(coeff);
                                indicesNeg = {i, j, k};
                            }
                        }
                        else if (numNeg == 2)
                        {
                            using std::abs;
                            auto coeffSorted = coeff;
                            std::sort(coeffSorted.begin(), coeffSorted.end());
                            const auto thisCost = coeffSorted[0]*coeffSorted[0] + coeffSorted[1]*coeffSorted[1];
                            foundValidNeg2Solution = true;

                            if (thisCost < costNeg2)
                            {
                                costNeg2 = thisCost;
                                bestCoefficientsNeg2 = std::move(coeff);
                                indicesNeg2 = {i, j, k};
                            }
                        }
                    }
                }
            }
        }

        if (foundValidPositiveSolution)
            return { std::move(indices), std::move(bestCoefficients), 3 };
        else if (foundValidNeg1Solution)
            return { std::move(indicesNeg), std::move(bestCoefficientsNeg), 3 };
        else if (foundValidNeg2Solution)
            return { std::move(indicesNeg2), std::move(bestCoefficientsNeg2), 3 };
        else
        {
            std::cout << "Found no decomposition of vector " << vec << " with\n";
            for (auto&& bb : b)
                std::cout << bb << " ;;; ";
            std::cout << std::endl;
            return { std::move(indices), std::move(bestCoefficients), 0 };
        }
    }
    else
        DUNE_THROW(Dune::NotImplemented, "Vector decomposition for dim = " << dim);
}

} // end namespace Dumux

#endif
