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
                const int maxP_ = 4096;
                assert(x1 != x2);
                
                x1_ = x1;
                x2_ = x2;
                
                if (p > 0) {
                    calcParameters_(x1, x2, y1, y2, m1, m2, p);
                    return;
                }
                
                // try to auto-detect a good value for p. we do this
                // by calculating the average second derivative of the
                // spline.
                p = 1;
                while (p <= maxP_)
                {
                    // calculate parameters for the current p
                    calcParameters_(x1, x2, y1, y2, m1, m2, p);
                    
                    if (isOk_(x1,x2,y1,y2,m1,m2))
                        break;

                    // try a bigger p next time
                    p = p*2;
                }

                if (p > 2 && p <= maxP_) {
                    int pLow = p/2;
                    int pHigh = p;
                    while (true) {
                        // calculate parameters for the current p
                        calcParameters_(x1, x2, y1, y2, m1, m2, p);
                        
                        if (isOk_(x1,x2,y1,y2,m1,m2))
                            pHigh = p;
                        else
                            pLow = p;
                        
                        if (pHigh - pLow < 4) {
                            p = pHigh;
                            calcParameters_(x1, x2, y1, y2, m1, m2, p);
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
            { return x1_ <= x && x <= x2_; };

        /*!
         * \brief Evaluate the polynomial at a given position.
         */
        Scalar eval(Scalar x) const
            {
                return a_ + 
                    b_*(x - x1_) +
                       c_*exp(p_*(x - x1_)) + 
                       d_*exp(-p_*(x - x1_));
            }

        /*!
         * \brief Evaluate the polynomial's derivative at a given position.
         */
        Scalar evalDerivative(Scalar x) const
            {
                return b_ + 
                    p_*( c_*exp( p_*(x - x1_)) - 
                         d_*exp(-p_*(x - x1_)) );
            }

        Scalar p() { return p_; }
        
    private:
        void calcParameters_(Scalar x1,
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
                
                a_ = dx*(beta + w*delta);
                b_ = alpha - beta + w*(gamma - delta);
                c_ = 0.5*(gamma - delta*exp(-p*dx))/(z - p);
                d_ = 0.5*(delta*exp(p*dx) - gamma)/(z - p);
                p_ = p;
            };

        bool isOk_(Scalar x1, Scalar x2,
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
                    ok = ok && isUnder_(x1, y1, m1/f);
                else 
                    ok = ok && isOver_(x1, y1, m1*f);

                if (m2 > mAvg)
                    ok = ok && isOver_(x2, y2, m2/f);
                else 
                    ok = ok && isUnder_(x2, y2, m2*f);

                return ok;
            }

        // return true iff the curve is unter the line which passes
        // through (x0, y0) with a slope of m
        bool isUnder_(Scalar x0, Scalar y0, Scalar m)
            {
                Scalar b = y0 - m*x0;

                Scalar x = x1_;
                Scalar inc = (x2_ - x1_)/(2 - 0.01);
                for (x += inc; x < x2_; x += inc) { 
                    // evaluate spline at position x
                    Scalar tmp = eval(x);

/*                    std::cerr << boost::format("isUnder_: %f > %f !\n")%
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
        bool isOver_(Scalar x0, Scalar y0, Scalar m)
            {
                Scalar b = y0 - m*x0;
                
                Scalar x = x1_;
                Scalar inc = (x2_ - x1_)/(2 - 0.01);
                for (x += inc; x < x2_; x += inc) { 
                    // evaluate spline at position x
                    Scalar tmp = eval(x);
                    
/*                    std::cerr << boost::format("isOver_: %f < %f !\n")%
                        (m*x + b)%
                        (tmp + fabs(tmp)*0.01);
*/

                    if (m*x + b > tmp - fabs(tmp)*0.01)
                        return false;
                }
                return true;
            };


        Scalar a_;
        Scalar b_;
        Scalar c_;
        Scalar d_;
        int    p_;

        Scalar x1_; 
        Scalar x2_; 
    };
}

#endif
