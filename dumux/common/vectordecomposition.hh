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

#ifndef DUMUX_COMMON_VECTORDECOMPOSITION_HH
#define DUMUX_COMMON_VECTORDECOMPOSITION_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include "math.hh"

namespace Dumux {

class VectorDecomposition
{

public:
    template<class DV>
    static std::tuple<std::vector<std::size_t>, std::vector<typename DV::value_type>, bool> calculateVectorDecomposition(const DV x, const std::vector<DV>& v)
    {
        std::size_t numVectors = v.size();
        static constexpr auto dim = DV::dimension;
        if(numVectors == 0)
            DUNE_THROW(Dune::InvalidStateException, "Can't perform decomposition without any given vectors!");

        if(isAligned(x, v[0]))
            return {std::vector<std::size_t>({0}), calculateCoefficients(x, v[0]), true};

        using Scalar = typename DV::value_type;
        std::vector<Scalar> coefficients(dim, std::numeric_limits<Scalar>::max());
        std::vector<std::size_t> indices;
        bool foundValidSolution = false;

        if constexpr (dim == 3)
        {
            for (std::size_t i = 0; i < numVectors - 2; ++i)
            {
                for (std::size_t j = i+1; j < numVectors - 1; ++j)
                {
                    for (std::size_t k = j+1; k < numVectors; ++k)
                    {
                        auto [coeff, valid] =  calculateCoefficients(x, v[i], v[j], v[k]);
                        using std::max_element; using std::move;
                        if(valid && (*max_element(coeff.begin(), coeff.end()) <  *max_element(coefficients.begin(), coefficients.end())))
                        {
                            coefficients = move(coeff);
                            indices = {i, j, k};
                            foundValidSolution = true;
                        }
                    }
                }
            }
        }
        else if constexpr (dim == 2)
        {
            for (std::size_t i = 0; i < numVectors - 1; ++i)
            {
                for (std::size_t j = i+1; j < numVectors; ++j)
                {
                    auto [coeff, valid] =  calculateCoefficients(x, v[i], v[j]);
                    using std::max_element; using std::move;
                    if(valid && (*max_element(coeff.begin(), coeff.end()) <  *max_element(coefficients.begin(), coefficients.end())))
                        {
                            coefficients = move(coeff);
                            indices = {i, j};
                            foundValidSolution = true;
                        }
                }
            }
        }

        return {indices, coefficients, foundValidSolution};
    }

private:
    //Cramer's rule for the solution of 2D system of equation
    template<class DV>
    static std::pair<std::vector<typename DV::value_type>, bool> calculateCoefficients(DV x, DV v1, DV v2)
    {
        std::vector<typename DV::value_type> coeff;
        const auto v1_norm = v1.two_norm();
        const auto v2_norm = v2.two_norm();

        v1 /= v1_norm;
        v2 /= v2_norm;
        const auto xNorm = x.two_norm();

        x /= xNorm;

        const auto A_f = v1[0]*v2[1] - v1[1]*v2[0];
        const auto A_f_1 = x[0]*v2[1] - x[1]*v2[0];
        const auto A_f_2 = v1[0]*x[1] - v1[1]*x[0];

        if(std::abs(A_f) < 1.0e-8)
            return {coeff, false};

        coeff.resize(2);
        coeff[0] = A_f_1 / A_f ;
        coeff[1] = A_f_2 / A_f ;

        for(auto& c : coeff)
            if(std::abs(c) < 1.0e-12)
                c = 0.0;

        coeff[0] *= xNorm/v1_norm;
        coeff[1] *= xNorm/v2_norm;

        if(coeff[0] >= -1.0e-30 && coeff[1] >= -1.0e-30)
            return {coeff, true};

        return {coeff, false};
    }

    template<class DV>
    static std::pair<std::vector<typename DV::value_type>, bool> calculateCoefficients(DV x, DV v1, DV v2, DV v3)
    {
        std::vector<typename DV::value_type> coeff;
        const auto v1_norm = v1.two_norm();
        const auto v2_norm = v2.two_norm();
        const auto v3_norm = v3.two_norm();
        const auto xNorm = x.two_norm();

        v1 /= v1_norm;
        v2 /= v2_norm;
        v3 /= v3_norm;

        x /= xNorm;

        const auto A_f = v1*crossProduct(v2,v3);
        const auto A_f_1 = x * crossProduct(v2,v3);
        const auto A_f_2 = v1 * crossProduct(x,v3);
        const auto A_f_3 = v1 * crossProduct(v2,x);

        if(std::abs(A_f) < 1.0e-8)
            return {coeff, false};

        coeff.resize(3);
        coeff[0] = A_f_1 / A_f ;
        coeff[1] = A_f_2 / A_f ;
        coeff[2] = A_f_3 / A_f ;

        for(auto& c : coeff)
            if(std::abs(c) < 1.0e-12)
                c = 0.0;

        coeff[0] *= xNorm/v1_norm;
        coeff[1] *= xNorm/v2_norm;
        coeff[2] *= xNorm/v3_norm;

        if(coeff[0] >=-1.0e-30 &&  coeff[1] >=-1.0e-30 && coeff[2] >=-1.0e-30)
            return {coeff, true};

        return {coeff, false};
    }

    template<class DV>
    static bool isAligned(DV x, DV v1)
    {
        const auto v1_norm = v1.two_norm();

        v1 /= v1_norm;
        const auto xNorm = x.two_norm();

        x /= xNorm;

        if(std::abs(v1*x - 1.0) < 1.0e-10)
            return true;

        return false;
    }

    template<class DV>
    static std::vector<typename DV::value_type> calculateCoefficients(const DV& x, const DV& v1)
    {
        std::vector<typename DV::value_type> coeff;
        coeff.resize(1);
        coeff[0] = std::abs(x*v1);
        coeff[0] /= v1.two_norm2();

        return coeff;
    }
};

} // end namespace Dumux

#endif
