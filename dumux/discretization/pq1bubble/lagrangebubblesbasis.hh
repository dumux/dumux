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
 * \ingroup PQ1BubbleDiscretization
 * \brief Evaluate basis with bubble function
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_BASIS_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_BASIS_HH

#include <array>
#include <numeric>

namespace Dumux {

class BubbleLocalBasis
{

public:
    //! \brief Evaluate all shape functions with bubble function
    template<class Geometry, class LocalBasis>
    static void evaluateFunctionWithBubble(const Geometry& geom,
                                           const LocalBasis& localBasis,
                                           const typename LocalBasis::Traits::DomainType& x,
                                           std::vector<typename LocalBasis::Traits::RangeType>& out)
    {
        localBasis.evaluateFunction(x, out);

        using ShapeValue = typename LocalBasis::Traits::RangeType;
        std::vector<ShapeValue> shapeValuesCenter;
        localBasis.evaluateFunction(geom.local(geom.center()), shapeValuesCenter);

        // We add the contribution of the bubble function
        ShapeValue bubbleFktValue(1.0);
        for(auto v : out)
            bubbleFktValue *= v;

        ShapeValue bubbleFktValueCenter(1.0);
        for(auto v : shapeValuesCenter)
            bubbleFktValueCenter *= v;

        for(int i=0; i<out.size(); i++)
            out[i] -= shapeValuesCenter[i]*bubbleFktValue/bubbleFktValueCenter;

        out.push_back(bubbleFktValue/bubbleFktValueCenter);
    }

    //! \brief Evaluate all shape functions with bubble function
    template<class Geometry, class LocalBasis>
    static void evaluateFunctionAndJacobianWithBubble(const Geometry& geom,
                                                      const LocalBasis& localBasis,
                                                      const typename LocalBasis::Traits::DomainType& x,
                                                      std::vector<typename LocalBasis::Traits::RangeType>& outF,
                                                      std::vector<typename LocalBasis::Traits::JacobianType>& outJ)
    {
        localBasis.evaluateFunction(x, outF);
        localBasis.evaluateJacobian(x, outJ);

        using ShapeValue = typename LocalBasis::Traits::RangeType;
        using ShapeJacobian = typename LocalBasis::Traits::JacobianType;

        std::vector<ShapeValue> shapeValuesCenter;
        localBasis.evaluateFunction(geom.local(geom.center()), shapeValuesCenter);

        // We add the contribution of the bubble function
        ShapeValue bubbleFktValue(1.0);
        for(auto v : outF)
            bubbleFktValue *= v;

        ShapeValue bubbleFktValueCenter(1.0);
        for(auto v : shapeValuesCenter)
            bubbleFktValueCenter *= v;

        ShapeJacobian bubbleJacobian(0.0);
        for(int i=0; i<outJ.size(); i++)
        {
            ShapeJacobian jacVal = outJ[i];
            for(int j=0; j<outF.size(); j++)
                if( i!= j)
                    jacVal *= outF[j];
            bubbleJacobian += jacVal;
        }
        bubbleJacobian /= bubbleFktValueCenter;

        for(int i=0; i<outF.size(); i++)
            outF[i] -= shapeValuesCenter[i]*bubbleFktValue/bubbleFktValueCenter;

        for(int i=0; i<outJ.size(); i++)
        {
            ShapeJacobian val = bubbleJacobian;
            val *= shapeValuesCenter[i];
            outJ[i] -= val;
        }

        outF.push_back(bubbleFktValue/bubbleFktValueCenter);
        outJ.push_back(bubbleJacobian);
    }
};

} // end namespace Dumux

#endif
