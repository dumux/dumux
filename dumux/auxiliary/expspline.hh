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
 * \brief An exponential spline.
 *
 * It has the same properties as a polynomial spline, but is more
 * robust. E.g.  it is much more likely that the slope of the curve is
 * monotonous.
 */
#ifndef DUNE_EXP_SPLINE_HH
#define DUNE_EXP_SPLINE_HH

#include <iostream>
#include <boost/format.hpp>

namespace Dune
{
    /*
     * \brief An exponential spline.
     *
     * It has the same properties as a polynomial spline, but is more
     * robust. E.g.  it is much more likely that the slope of the curve is
     * monotonous.
     */
    template<class ScalarT>
    class ExpSpline
    {
        typedef ScalarT Scalar;

    public:
        ExpSpline(Scalar x1, 
                  Scalar x2, 
                  Scalar y1, 
                  Scalar y2, 
                  Scalar m1, 
                  Scalar m2,
                  int p = 0)
            {
                const int _maxP = 4096;
                assert(x1 != x2);
                
                _x1 = x1;
                _x2 = x2;
                
                if (p > 0) {
                    _calcParameters(x1, x2, y1, y2, m1, m2, p);
                    return;
                }
                
                // try to auto-detect a good value for p. we do this
                // by calculating the average second derivative of the
                // spline.
                p = 1;
                while (p <= _maxP)
                {
                    // calculate parameters for the current p
                    _calcParameters(x1, x2, y1, y2, m1, m2, p);
                    
                    if (_isOk(x1,x2,y1,y2,m1,m2))
                        break;

                    // try a bigger p next time
                    p = p*2;
                }

                if (p > 2 && p <= _maxP) {
                    int pLow = p/2;
                    int pHigh = p;
                    while (true) {
                        // calculate parameters for the current p
                        _calcParameters(x1, x2, y1, y2, m1, m2, p);
                        
                        if (_isOk(x1,x2,y1,y2,m1,m2))
                            pHigh = p;
                        else
                            pLow = p;
                        
                        if (pHigh - pLow < 4) {
                            p = pHigh;
                            _calcParameters(x1, x2, y1, y2, m1, m2, p);
                            break;
                        }
                        
                        // use the middle p for the next round
                        p = (pHigh + pLow)/2;
                    }
                }
                
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
                return _a + 
                    _b*(x - _x1) +
                       _c*exp(_p*(x - _x1)) + 
                       _d*exp(-_p*(x - _x1));
            }

        /*!
         * \brief Evaluate the polynomial's derivative at a given position.
         */
        Scalar evalDerivative(Scalar x) const
            {
                return _b + 
                    _p*( _c*exp( _p*(x - _x1)) - 
                         _d*exp(-_p*(x - _x1)) );
            }

        Scalar p() { return _p; }
        
    private:
        void _calcParameters(Scalar x1,
                             Scalar x2,
                             Scalar y1,
                             Scalar y2,
                             Scalar m1,
                             Scalar m2,
                             int p)
            {
                // calculate the coefficents of the exponential spline.
                Scalar dx = x2 - x1;
                Scalar dy = y2 - y1;
                
                Scalar z = sinh(p*dx)/dx;
                Scalar w = -z/(z - p);
                
                Scalar alpha = y2/dx;
                Scalar beta = y1/dx;
                
                Scalar v = (p*cosh(p*dx) - z)/(z - p);
                
                Scalar gamma = (m1 + v*m2 - (v + 1)*dy/dx)/(v*v - 1);
                Scalar delta = (-v*m1 - m2 + (v + 1)*dy/dx)/(v*v - 1);
                
                _a = dx*(beta + w*delta);
                _b = alpha - beta + w*(gamma - delta);
                _c = 0.5*(gamma - delta*exp(-p*dx))/(z - p);
                _d = 0.5*(delta*exp(p*dx) - gamma)/(z - p);
                _p = p;
            };

        bool _isOk(Scalar x1, Scalar x2,
                   Scalar y1, Scalar y2,
                   Scalar m1, Scalar m2)
            {
                if (x1 > x2) {
                    std::swap(x1,x2);
                    std::swap(y1,y2);
                    std::swap(m1,m2);
                }
   
                Scalar mAvg = (y2 - y1)/(x2 - x1);

                bool ok = true;
                Scalar f = 1.01;

                if (m1 > mAvg)
                    ok = ok && _isUnder(x1, y1, m1/f);
                else 
                    ok = ok && _isOver(x1, y1, m1*f);

                if (m2 > mAvg)
                    ok = ok && _isOver(x2, y2, m2/f);
                else 
                    ok = ok && _isUnder(x2, y2, m2*f);

                return ok;
            }

        // return true iff the curve is unter the line which passes
        // through (x0, y0) with a slope of m
        bool _isUnder(Scalar x0, Scalar y0, Scalar m)
            {
                Scalar b = y0 - m*x0;

                Scalar x = _x1;
                Scalar inc = (_x2 - _x1)/(2 - 0.01);
                for (x += inc; x < _x2; x += inc) { 
                    // evaluate spline at position x
                    Scalar tmp = eval(x);

/*                    std::cerr << boost::format("_isUnder: %f > %f !\n")%
                        (m*x + b)%
                        (tmp + fabs(tmp)*0.01);
*/
                    if (m*x + b < tmp + fabs(tmp)*0.01)
                        return false;
                }
                return true;
            };

        // return true iff the curve is unter the line which passes
        // through (x0, y0) with a slope of m
        bool _isOver(Scalar x0, Scalar y0, Scalar m)
            {
                Scalar b = y0 - m*x0;
                
                Scalar x = _x1;
                Scalar inc = (_x2 - _x1)/(2 - 0.01);
                for (x += inc; x < _x2; x += inc) { 
                    // evaluate spline at position x
                    Scalar tmp = eval(x);
                    
/*                    std::cerr << boost::format("_isOver: %f < %f !\n")%
                        (m*x + b)%
                        (tmp + fabs(tmp)*0.01);
*/

                    if (m*x + b > tmp - fabs(tmp)*0.01)
                        return false;
                }
                return true;
            };


        Scalar _a;
        Scalar _b;
        Scalar _c;
        Scalar _d;
        int    _p;

        Scalar _x1; 
        Scalar _x2; 
    };
}

#endif
