/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief A 3rd order polynomial p(x) for which, given two points x1 and x2!=x1, 
 *        the following hold:
 *        p(x1) = y1
 *        p(x2) = y2
 *        p'(x1) = m1
 *        p'(x2) = m2
 *
 * or any given y1, y2, m1, m2.
 */
#ifndef DUNE_SPLINE_HH
#define DUNE_SPLINE_HH

namespace Dune
{
    /*
     * \brief A 2rd order polynomial p(x) = a x^3 + b x^2 + c x + d for 
     *        which, given two distinct points x1 and x2, the following holds:
     *        p(x1) = y1
     *        p(x2) = y2
     *        p'(x1) = m1
     *        p'(x2) = m2
     *
     * for any given y1, y2, m1, m2. This class is quite useful for
     * material laws such as Parker-Lenhard to avoid discontinuities
     * in the first derivative which result in numerical problems for
     * the newton solver.
     */
    template<class ScalarT>
    class Spline
    {
        typedef ScalarT Scalar;

    public:
        Spline(Scalar x1, 
               Scalar x2, 
               Scalar y1, 
               Scalar y2, 
               Scalar m1, 
               Scalar m2)
            {
                assert(x1 != x2);
                
                _x1 = x1;
                _x2 = x2;
                
                // calculate the coefficents of the polynomial. this
                // is pretty cumbersome, since we basically solve a
                // 4x4 matrix analytically here.
                Scalar tmpLeft  = 2*(x1*x1*x1 - x2*x2*x2) - 3*(x1 - x2)*(x1*x1 + x2*x2);
                Scalar tmpRight = 2*(y1 - y2) - (x1 - x2)*(m1 + m2);
                _a = tmpRight/tmpLeft;
                
                tmpLeft = 2*(x1 - x2);
                tmpRight = m1 - m2 - (3*(x1*x1 - x2*x2)*_a);
                _b = tmpRight/tmpLeft;

                tmpRight = m1 - (3*x1*x1*_a + 2*x1*_b);
                tmpLeft = 1;
                _c = tmpRight/tmpLeft;
                
                tmpRight = y1 - (x1*x1*x1*_a + x1*x1*_b + x1*_c);
                tmpLeft  = 1;
                _d = tmpRight/tmpLeft;

                assert(fabs(eval(x1) - y1) < 1e-5);
                assert(fabs(eval(x2) - y2) < 1e-5);
                assert(fabs(evalDerivative(x1) - m1) < 1e-5);
                assert(fabs(evalDerivative(x2) - m2) < 1e-5);
            }
        
        /*!
         * \brief Return true iff the given x is in range [x1, x2].
         */
        bool applies(Scalar x) const
            { return _x1 <= x && x <= _x2; };

        /*!
         * \brief Evaluate the polynomial at a given position.
         */
        Scalar eval(Scalar x) const
            {
                return ((_a*x + _b)*x + _c)*x + _d;
            }

        /*!
         * \brief Evaluate the polynomial's derivative at a given position.
         */
        Scalar evalDerivative(Scalar x) const
            {
                return (3*_a*x + 2*_b)*x + _c;
            }
        
    private:
        Scalar _a;
        Scalar _b;
        Scalar _c;
        Scalar _d;

        Scalar _x1; 
        Scalar _x2; 
    };
}

#endif
