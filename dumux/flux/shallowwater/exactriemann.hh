// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup ShallowWater
 * \brief Function to compute the Riemann flux at the interface
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_EXACT_RIEMANN_HH
#define DUMUX_FLUX_SHALLOW_WATER_EXACT_RIEMANN_HH

#include <array>
#include <cmath>

namespace Dumux {
namespace ShallowWater{

/*!
 * \ingroup ShallowWater
 * \brief Exact Riemann solver for Shallow water equations.
 *
 * This Riemann solver is described in the book
 * "Shock-capturing methods for free-surface shallow flows"
 * from Toro, 2001. We keep the notation for the variables
 * after Toro.
 *
 * \param dl water depth on the left side
 * \param dr water depth on the right side
 * \param ul veloctiyX on the left side
 * \param ur velocityX on the right side
 * \param vl velocityY on the left side
 * \param vr velocityY on the right side
 * \param grav gravity constant
 */
template<class Scalar>
std::array<Scalar,3> exactRiemann(const Scalar& dl,
                                  const Scalar& dr,
                                  const Scalar& ul,
                                  const Scalar& ur,
                                  const Scalar& vl,
                                  const Scalar& vr,
                                  const Scalar& grav)
{

        Scalar fl = 0.0; //function fl
        Scalar fr = 0.0; //function fr
        Scalar fld = 0.0; //function fl derivative
        Scalar frd = 0.0; //function fr derivative
        Scalar us = 0.0; //ustar
        Scalar cl = 0.0; //celerity sqrt(g*h) left side
        Scalar cr = 0.0; //celerity sqrt(g*h) right side
        Scalar c = 0.0; //celerity
        Scalar cs = 0.0; //celerity
        Scalar s = 0.0; //sample point
        Scalar dcrit = 0.0; //critical water depth
        Scalar dmin = 0.0; //minimum depth
        Scalar gel = 0.0; //gravity component
        Scalar ger = 0.0; //gravity component
        Scalar ges = 0.0; //gravity component
        Scalar ds = 0.0; //water depth
        Scalar d0 = 0.0; //start point for Newton method
        Scalar cha = 0.0; //variable for Newton method
        Scalar sl = 0.0; //shock left
        Scalar sr = 0.0; //shock right
        Scalar shl = 0.0; //wave speed ad head left
        Scalar shr = 0.0; //wave speed ad head right
        Scalar stl = 0.0; //wave speed at tail left
        Scalar str = 0.0; // wave speed at tail right
        int maxsteps = 200; //maximum steps for the Newton method
        Scalar tol = 1.0E-12; //tolerance for the Newton method
        Scalar d = 0.0; // water depth
        Scalar u = 0.0; //velocityX
        Scalar v = 0.0; //celocityY
        Scalar ql = 0.0; //flux left
        Scalar qr = 0.0; //flux right
        Scalar ssl = 0.0; //wave speed ssl
        Scalar ssr = 0.0; //wave speed ssr

        using std::sqrt;
        using std::max;
        using std::min;
        using std::pow;
        using std::abs;

        //================== Exact Riemann solver
        cl = sqrt(grav * max(dl,0.0));
        cr = sqrt(grav * max(dr,0.0));
        dcrit = (ur - ul) - 2.0*(cl+cr);

        //we use s as 0.0 since we are at the x = 0
        s = 0.0;

        // dry case
        if ((dl <= 0)||(dr <= 0)||(dcrit >= 0))
        {
            //left side is dry
            if(dl <= 0)
            {
                shr = ur + cr;

                if (s >= shr)
                {
                    d = dr;
                    u = ur;
                    v = vr;

                }

                else
                {

                    str = ur -2.0*cr;

                    if(s >= str){
                        u = (ur - 2.0 *cr + 2.0*s)/3.0;
                        c = (-ur + 2.0*cr +s )/3.0;
                        d = c*c/grav;
                        v = vr;

                    }

                    else
                    {
                        d = dl;
                        u = ul;
                        v = vl;
                    }
                }
            }

            else
            {
                //right side is dry
                if(dr <= 0)
                {
                    shl = ul -cl;

                    if(s <= shl)
                    {
                        d = dl;
                        u = ul;
                        v = vl;
                    }

                    else
                    {
                        stl = ul + 2.0 *cl;

                        if ( s <= stl)
                        {
                            u = (ul + 2.0 *cl + 2.0*s)/3.0;
                            c = (ul + 2.0*cl -s)/3.0;
                            d = c*c/grav;
                            v = vl;
                        }

                        else
                        {
                            d = dr;
                            u = ur;
                            v = vr;
                        }
                    }
                }

                //middle state is dry
                else
                {
                    shl = ul - cl;
                    ssl = ul + 2.0 * cl;
                    ssr = ur - 2.0 * cr;
                    shr = ur + cr;

                    if(s <= shl)
                    {
                        d = dl;
                        u = ul;
                        v = vl;
                    }

                    if((s > shl)&&(s <= ssl))
                    {
                        u = (ul + 2.0 * cl + 2.0*s)/3.0;
                        c = (ul + 2.0*cl-s)/3.0;
                        d = c*c/grav;
                        v = vl;
                    }

                    if((s > ssl)and (s<= ssr))
                    {
                        d = 0.0;
                        u = 0.0;
                        v = 0.0;
                    }

                    if ((s > ssr ) &&(s <= shr))
                    {
                        u = (ur - 2.0*cr + 2.0 *s)/3.0;
                        c = (-ur + 2.0*cr +s)/3.0;
                        d = c*c/grav;
                        v = vr;
                    }

                    if(s > shr)
                    {
                        d = dr;
                        u = ur;
                        v = vr;
                    }
                }
            }
        }

        // wet case
        else
        {
            //Get a starting value for the Newton-Raphson method
            ds = (1.0/grav) * pow(0.5 * (cl + cr) - 0.25 * (ul - ur),2.0);

            //default is Two-Rarefaction  as starting value, if ds > dmin we use two-shock
            dmin = min(dl,dr);

            if (ds > dmin)
            {
                gel = sqrt(0.5*grav * (ds+dl)/(ds*dl));
                ger = sqrt(0.5*grav * (ds+dr)/(ds*dr));
                ds = (gel * dl + ger * dr - (ur-ul))/(gel + ger);
            }

            //Start Newton-Raphson loop ds = hstar
            ds = max(ds,tol);
            d0 = ds;

            for (int i=0; i <= maxsteps; ++i)
            {
                //Compute fl,fr,fld,frd

                //first left side
                if (ds <= dl)
                {
                    c = sqrt(grav * ds);
                    fl = 2.0 * (c-cl);
                    fld = grav/c;
                }

                else
                {
                    ges = sqrt(0.5 * grav * (ds+dl)/(ds*dl));
                    fl = (ds - dl) * ges;
                    fld = ges - 0.25 * grav * (ds - dl)/(ges * ds * ds);
                }

                //second right side
                if (ds <= dr)
                {
                    c = sqrt(grav * ds);
                    fr = 2.0 * (c-cr);
                    frd = grav/c;
                }

                else
                {
                    ges = sqrt(0.5 * grav * (ds+dr)/(ds*dr));
                    fr = (ds - dr) * ges;
                    frd = ges - 0.25 * grav * (ds - dr)/(ges * ds * ds);
                }

                ds -= (fl + fr + ur - ul)/(fld + frd);
                cha =  abs(ds - d0)/(0.5*(ds+d0));

                if (cha <= tol)
                {
                    break;
                }

                if (ds < 0.0)
                {
                    ds = tol;
                }

                d0 = ds;
            }

            //compute ustar (us)
            if (ds <= dl)
            {
                c = sqrt(grav * ds);
                fl = 2.0 * (c-cl);
            }

            else
            {
                ges = sqrt(0.5 * grav * (ds+dl)/(ds*dl));
                fl = (ds - dl) * ges;
            }

            if (ds <= dr)
            {
                c = sqrt(grav * ds);
                fr = 2.0 * (c-cr);
            }

            else
            {
                ges = sqrt(0.5 * grav * (ds+dr)/(ds*dr));
                fr = (ds - dr) * ges;
            }

            //exact Riemann solver sstar = ustar
            us = 0.5 * (ul + ur) + 0.5 * (fr-fl);
            cs = sqrt(grav * ds);

            /***********************  computation of u and d *******************/
            //left wave
            if (s <= us)
            {
                v = vl;
                //left shock
                if (ds >= dl)
                {
                    ql = sqrt((ds + dl)*ds/(2.0*dl*dl));
                    sl = ul - cl * ql;

                    //left side of shock
                    if(s <= sl)
                    {
                        d = dl;
                        u = ul;
                    }

                    //right side of shock
                    else
                    {
                        d = ds;
                        u = us;
                    }
                }

                //left rarefaction
                else
                {
                    shl = ul -cl;

                    //right side of rarefaction
                    if (s <= shl)
                    {
                        d = dl;
                        u = ul;
                    }

                    else
                    {
                        stl = us -cs;

                        //inside the rarefaction
                        if(s <= stl)
                        {
                            u = (ul + 2.0*cl + 2.0 * s)/3.0;
                            c = (ul + 2.0* cl -s)/3.0;
                            d = c*c/grav;
                        }

                        //inside the star region
                        else
                        {
                            d = ds;
                            u = us;
                        }
                    }
                }
            }

            //right wave
            else
            {
                v = vr;

                //right shock
                if(ds >= dr)
                {
                    qr = sqrt((ds +dr)* ds /(2.0*dr*dr));
                    sr = ur  + cr * qr;

                    //right side of shock
                    if(s >= sr)
                    {
                        d = dr;
                        u = ur;
                    }

                    //left side of shock
                    else
                    {
                        d = ds;
                        u = us;
                    }

                //right rarefaction
                }

                else
                {
                    shr  = ur + cr;

                    //right side of Rarefaction
                    if (s >= shr)
                    {
                        d = dr;
                        u = ur;
                    }

                    else
                    {
                        str = us + cs;

                        //inside the rarefaction
                        if(s>=str)
                        {
                            u = (ur -2.0 * cr + 2.0 * s)/3.0;
                            c = (-ur + 2.0 * cr + s)/3.0;
                            d = c*c/grav;

                        }

                        //inside the star region
                        else
                        {
                            d = ds;
                            u = us;
                        }
                    }
                }
            }
        }
        //============================== Flux computation ===================
        std::array<Scalar,3> flux;

        flux[0] = d*u;
        flux[1] = d * u * u + 0.5 * grav * d * d;
        flux[2] = d*u*v;

        return flux;
}
} // end namespace ShallowWater
} // end namespace Dumux

#endif
