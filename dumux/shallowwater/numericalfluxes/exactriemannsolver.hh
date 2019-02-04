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
 * \ingroup SweModel
 * \brief Exact Riemann solver for shallow water equations using
 *        two-point flux approximation.
 */
#ifndef DUMUX_SHALLOWWATER_NUMERICALFLUXES_EXACTRIEMANNSOLVER_HH
#define DUMUX_SHALLOWWATER_NUMERICALFLUXES_EXACTRIEMANNSOLVER_HH


#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>


namespace Dumux
{

template<class Scalar>

    //We use the notation after Toro dl = hl,
    //Toro, E.F. 2001, "Shock-Capturing Methods for Free-Surface Shallow Flows", 1st edition, John Wiley & Sons.

    //typedef Scalar ;
    void computeExactRiemann(Scalar *flux, Scalar dl, Scalar dr, Scalar ul, Scalar ur,Scalar vl, Scalar vr,Scalar grav){

        Scalar fl = 0.0;
        Scalar fr = 0.0;
        Scalar fld = 0.0;
        Scalar frd = 0.0;
        Scalar us = 0.0;
        Scalar hstar = 0.0;
        Scalar cl = 0.0; //cl
        Scalar cr = 0.0; //Cr
        Scalar c = 0.0;
        Scalar cs = 0.0;
        Scalar f = 0.0;
        Scalar s = 0.0;
        Scalar fd = 0.0;
        Scalar dcrit = 0.0;
        Scalar dmin = 0.0;
        Scalar gel = 0.0;
        Scalar ger = 0.0;
        Scalar ges = 0.0;
        Scalar ds = 0.0;
        Scalar d0 = 0.0;
        Scalar cha = 0.0;
        Scalar sl = 0.0;
        Scalar sr = 0.0;
        Scalar shl = 0.0;
        Scalar shr = 0.0;
        Scalar stl = 0.0;
        Scalar str = 0.0;
        Scalar utmp = 0.0;
        int maxsteps = 200;
        Scalar tol = 1.0E-12;
        Scalar d = 0.0;
        Scalar u = 0.0;
        Scalar ql = 0.0;
        Scalar qr = 0.0;
        Scalar v = 0.0;
        Scalar ds_start;
        Scalar ds_start2;
        Scalar ssl = 0.0;
        Scalar ssr = 0.0;

        //================== Exact Riemann solver
        dl = std::max(dl,0.0);
        dr = std::max(dr,0.0);

        cl = std::sqrt(grav * dl);
        cr = std::sqrt(grav * dr);
        dcrit = (ur - ul) - 2.0*(cl+cr);

        //we use s as 0.0 since we are at the x = 0
        s = 0.0;

        //--- dry case ---/
        if ((dl <= 0)||(dr <= 0)||(dcrit >= 0)){

            //left side is dry
            if(dl <= 0){
                shr = ur + cr;

                if (s >= shr){
                    d = dr;
                    u = ur;
                    v = vr;
                }else{

                    str = ur -2.0*cr;

                    if(s >= str){
                        u = (ur - 2.0 *cr + 2.0*s)/3.0;
                        c = (-ur + 2.0*cr +s )/3.0;
                        d = c*c/grav;
                        v = vr;

                    }else{
                        d = dl;
                        u = ul;
                        v = vl;
                    }
                }
            }else{
                //right side is dry
                if(dr <= 0){
                    shl = ul -cl;

                    if(s <= shl){
                        d = dl;
                        u = ul;
                        v = vl;
                    }else{
                        stl = ul + 2.0 *cl;

                        if ( s <= stl){
                            u = (ul + 2.0 *cl + 2.0*s)/3.0;
                            c = (ul + 2.0*cl -s)/3.0;
                            d = c*c/grav;
                            v = vl;
                        }else{
                            d = dr;
                            u = ur;
                            v = vr;
                        }
                    }
                //middle state is dry
                }else{
                    shl = ul - cl;
                    ssl = ul + 2.0 * cl;
                    ssr = ur - 2.0 * cr;
                    shr = ur + cr;

                    if(s <= shl){
                        d = dl;
                        u = ul;
                        v = vl;
                    }
                    if((s > shl)&&(s <= ssl)){
                        u = (ul + 2.0 * cl + 2.0*s)/3.0;
                        c = (ul + 2.0*cl-s)/3.0;
                        d = c*c/grav;
                        v = vl;
                    }
                    if((s > ssl)and (s<= ssr)){
                        d = 0.0;
                        u = 0.0;
                        v = 0.0;
                    }
                    if ((s > ssr ) &&(s <= shr)){
                        u = (ur - 2.0*cr + 2.0 *s)/3.0;
                        c = (-ur + 2.0*cr +s)/3.0;
                        d = c*c/grav;
                        v = vr;
                    }
                    if(s > shr){
                        d = dr;
                        u = ur;
                        v = vr;
                    }
                }
            }
        }else{
            //---- wet case ---/

            //Get a starting value for the Newton-Raphson method
            ds = (1.0/grav) * std::pow(0.5 * (cl + cr) - 0.25 * (ul - ur),2.0);
            ds_start2 = ds;

            //default is Two-Rarefaction  as starting value, if ds > dmin we use two-shock
            dmin = std::min(dl,dr);
            if (ds > dmin){
                gel = std::sqrt(0.5*grav * (ds+dl)/(ds*dl));
                ger = std::sqrt(0.5*grav * (ds+dr)/(ds*dr));
                ds = (gel * dl + ger * dr - (ur-ul))/(gel + ger);
            }

            //Start Newton-Raphson loop ds = hstar
            ds = std::max(ds,tol);
            d0 = ds;
            ds_start = ds;

            for (int i=0; i <= maxsteps; ++i){

                //Compute fl,fr,fld,frd
                //-- first left side
                if (ds <= dl){
                    c = std::sqrt(grav * ds);
                    fl = 2.0 * (c-cl);
                    fld = grav/c;
                }else{
                    ges = std::sqrt(0.5 * grav * (ds+dl)/(ds*dl));
                    fl = (ds - dl) * ges;
                    fld = ges - 0.25 * grav * (ds - dl)/(ges * ds * ds);
                }
                //-- second right side
                if (ds <= dr){
                    c = std::sqrt(grav * ds);
                    fr = 2.0 * (c-cr);
                    frd = grav/c;
                }else{
                    ges = std::sqrt(0.5 * grav * (ds+dr)/(ds*dr));
                    fr = (ds - dr) * ges;
                    frd = ges - 0.25 * grav * (ds - dr)/(ges * ds * ds);
                }

                ds -= (fl + fr + ur - ul)/(fld + frd);
                cha =  std::abs(ds - d0)/(0.5*(ds+d0));

                if (cha <= tol){
                    break;
                }
                if (ds < 0.0){
                    ds = tol;
                }
                d0 = ds;
                /*if (i == maxsteps){
                    std::cout << "\nWarning exact Riemann solver Netwon-Raphson not converged" << std::endl;
                    std::cout << "ds start " << ds_start << " ds " << ds << std::endl;
                    std::cout << "ds start2 " << ds_start2 << " dmin " << dmin << std::endl;
                    std::cout << "dl  " << dl << " dr " << dr << std::endl;
                    std::cout << "dcrit " << dcrit << std::endl;
                }*/
            }

            //compute ustar (us)
            if (ds <= dl){
                c = std::sqrt(grav * ds);
                fl = 2.0 * (c-cl);
            }else{
                ges = std::sqrt(0.5 * grav * (ds+dl)/(ds*dl));
                fl = (ds - dl) * ges;
            }
            if (ds <= dr){
                c = std::sqrt(grav * ds);
                fr = 2.0 * (c-cr);
            }else{
                ges = std::sqrt(0.5 * grav * (ds+dr)/(ds*dr));
                fr = (ds - dr) * ges;
            }

            //exact Riemann solver sstar = ustar
            us = 0.5 * (ul + ur) + 0.5 * (fr-fl);
            cs = std::sqrt(grav * ds);

            /***********************  computation of u and d *******************/
            //left wave
            if (s <= us){
                v = vl;
                //left shock
                if (ds >= dl){
                    ql = std::sqrt((ds + dl)*ds/(2.0*dl*dl));
                    sl = ul - cl * ql;

                    //left side of shock
                    if(s <= sl){
                        d = dl;
                        u = ul;
                    //right side of shock
                    }else{
                        d = ds;
                        u = us;
                    }
                //left rarefaction
                }else{
                    shl = ul -cl;

                    //right side of rarefaction
                    if (s <= shl){
                        d = dl;
                        u = ul;
                    }else{
                        stl = us -cs;

                        //inside the rarefaction
                        if(s <= stl){
                            u = (ul + 2.0*cl + 2.0 * s)/3.0;
                            c = (ul + 2.0* cl -s)/3.0;
                            d = c*c/grav;

                        //inside the star region
                        }else{
                            d = ds;
                            u = us;
                        }
                    }
                }
            //right wave
            }else{
                //std::cout << "right wave" << std::endl;
                v = vr;
                //right shock
                if(ds >= dr){
                    qr = std::sqrt((ds +dr)* ds /(2.0*dr*dr));
                    sr = ur  + cr * qr;

                    //std::cout << "sr " << sr << std::endl;
                    //right side of shock
                    if(s >= sr){
                        d = dr;
                        u = ur;
                    //left side of shock
                    }else{
                        d = ds;
                        u = us;
                    }

                //right rarefaction
                }else{
                    shr  = ur + cr;
                    //std::cout << "shr " << shr << " s " << s << std::endl;

                    //right side of Rarefaction
                    if (s >= shr){
                        //std::cout << "debug " << std::endl;
                        d = dr;
                        u = ur;
                    }else{
                        str = us + cs;
                        //inside the rarefaction
                        if(s>=str){
                            u = (ur -2.0 * cr + 2.0 * s)/3.0;
                            c = (-ur + 2.0 * cr + s)/3.0;
                            d = c*c/grav;

                        //inside the star region
                        }else{
                            d = ds;
                            u = us;
                        }
                    }
                }
            }
        }
        //============================== Flux computation ===================0

        //now we have d and u
        //std::cout << "d " << d <<  " u " << u << std::endl;
        flux[0] = d*u;
        flux[1] = d * u * u + 0.5 * grav * d * d;
        flux[2] = d*u*v;
}
} // end namespace Dumux

#endif
