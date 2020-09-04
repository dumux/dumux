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
 * \ingroup ShallowWaterFlux
 * \brief Function to compute the Riemann flux at the interface
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_EXACT_RIEMANN_HH
#define DUMUX_FLUX_SHALLOW_WATER_EXACT_RIEMANN_HH

#include <array>
#include <cmath>

namespace Dumux {
namespace ShallowWater{

template<typename Scalar>
struct RiemannSolution {
    std::array<Scalar,3> flux;
    Scalar waterDepth;
    Scalar velocityX;
    Scalar velocityY;
};


/*!
 * \ingroup ShallowWaterFlux
 * \brief Exact Riemann solver for the shallow water equations.
 *
 * The flux of the 2D shallow water equations must be rotated
 * to a 1D problem before the Riemann solver can be applied.
 * The computed water flux is given in m^2/s, the momentum
 * fluxes are given in m^3/s^2.
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
 * \param s sample point (default = 0 since x = 0 for flux computation)
 */
template<class Scalar>
RiemannSolution<Scalar> exactRiemann(const Scalar dl,
                                     const Scalar dr,
                                     const Scalar ul,
                                     const Scalar ur,
                                     const Scalar vl,
                                     const Scalar vr,
                                     const Scalar grav,
                                     const Scalar s = 0.0)
{
    RiemannSolution<Scalar> sol;
    sol.waterDepth = 0.0; // d (Toro)
    sol.velocityX = 0.0; // u (Toro)
    sol.velocityY = 0.0; // v (Toro)

    constexpr int maxsteps = 200; //maximum steps for the Newton method
    constexpr Scalar tol = 1.0E-12; //tolerance for the Newton method

    using std::sqrt;
    using std::max;
    using std::min;
    using std::abs;

    //================== Exact Riemann solver
    const auto cl = sqrt(grav * max(dl, 0.0)); //celerity sqrt(g*h) left side
    const auto cr = sqrt(grav * max(dr, 0.0)); //celerity sqrt(g*h) right side
    const auto dcrit = (ur - ul) - 2.0*(cl + cr); //critical water depth

    // dry case
    if ((dl <= 0) || (dr <= 0) || (dcrit >= 0))
    {
        //left side is dry
        if(dl <= 0)
        {
            const auto shr = ur + cr; //wave speed at head right

            if (s >= shr)
            {
                sol.waterDepth = dr;
                sol.velocityX = ur;
                sol.velocityY = vr;

            }

            else
            {

                const auto str = ur - 2.0*cr; //wave speed at tail right

                if(s >= str){
                    sol.velocityX = (ur - 2.0 *cr + 2.0*s)/3.0;
                    const auto c = (-ur + 2.0*cr +s )/3.0;
                    sol.waterDepth = c*c/grav;
                    sol.velocityY = vr;

                }

                else
                {
                    sol.waterDepth = dl;
                    sol.velocityX = ul;
                    sol.velocityY = vl;
                }
            }
        }

        else
        {
            //right side is dry
            if(dr <= 0)
            {
                const auto shl = ul-cl; //wave speed at head left

                if(s <= shl)
                {
                    sol.waterDepth = dl;
                    sol.velocityX = ul;
                    sol.velocityY = vl;
                }

                else
                {
                    const auto stl = ul + 2.0 *cl; //wave speed at tail left

                    if ( s <= stl)
                    {
                        sol.velocityX = (ul + 2.0 *cl + 2.0*s)/3.0;
                        const auto c = (ul + 2.0*cl - s)/3.0;
                        sol.waterDepth = c*c/grav;
                        sol.velocityY = vl;
                    }

                    else
                    {
                        sol.waterDepth = dr;
                        sol.velocityX = ur;
                        sol.velocityY = vr;
                    }
                }
            }

            //middle state is dry
            else
            {
                const auto shl = ul - cl; //wave speed at head left
                const auto shr = ur + cr; //wave speed at head right
                const auto ssl = ul + 2.0 * cl;
                const auto ssr = ur - 2.0 * cr;

                if(s <= shl)
                {
                    sol.waterDepth = dl;
                    sol.velocityX = ul;
                    sol.velocityY = vl;
                }

                if((s > shl) && (s <= ssl))
                {
                    sol.velocityX = (ul + 2.0 * cl + 2.0*s)/3.0;
                    const auto c = (ul + 2.0*cl - s)/3.0;
                    sol.waterDepth = c*c/grav;
                    sol.velocityY = vl;
                }

                if((s > ssl) && (s <= ssr))
                {
                    sol.waterDepth = 0.0;
                    sol.velocityX = 0.0;
                    sol.velocityY = 0.0;
                }

                if ((s > ssr ) && (s <= shr))
                {
                    sol.velocityX = (ur - 2.0*cr + 2.0*s)/3.0;
                    const auto c = (-ur + 2.0*cr + s)/3.0;
                    sol.waterDepth = c*c/grav;
                    sol.velocityY = vr;
                }

                if(s > shr)
                {
                    sol.waterDepth = dr;
                    sol.velocityX = ur;
                    sol.velocityY = vr;
                }
            }
        }
    }

    // wet case
    else
    {
        // Get a starting value for the Newton-Raphson method
        const auto tmp = 0.5 * (cl + cr) - 0.25 * (ul - ur);
        auto ds = (1.0/grav) * tmp * tmp; // water depth

        //default is Two-Rarefaction  as starting value, if ds > dmin we use two-shock
        const auto dmin = min(dl, dr); // minimum depth

        if (ds > dmin)
        {
            const auto gel = sqrt(0.5*grav * (ds+dl)/(ds*dl));
            const auto ger = sqrt(0.5*grav * (ds+dr)/(ds*dr));
            ds = (gel * dl + ger * dr - (ur-ul))/(gel + ger);
        }

        //Start Newton-Raphson loop ds = hstar
        ds = max(ds, tol);
        auto d0 = ds;

        Scalar fl = 0.0, fr = 0.0, fld = 0.0, frd = 0.0;

        for (int i=0; i <= maxsteps; ++i)
        {
            //Compute fl,fr,fld,frd

            //first left side
            if (ds <= dl)
            {
                const auto c = sqrt(grav * ds);
                fl = 2.0 * (c-cl);
                fld = grav/c;
            }

            else
            {
                const auto ges = sqrt(0.5 * grav * (ds+dl)/(ds*dl));
                fl = (ds - dl) * ges;
                fld = ges - 0.25 * grav * (ds - dl)/(ges * ds * ds);
            }

            //second right side
            if (ds <= dr)
            {
                const auto c = sqrt(grav * ds);
                fr = 2.0 * (c-cr);
                frd = grav/c;
            }

            else
            {
                const auto ges = sqrt(0.5 * grav * (ds+dr)/(ds*dr));
                fr = (ds - dr) * ges;
                frd = ges - 0.25 * grav * (ds - dr)/(ges * ds * ds);
            }

            ds -= (fl + fr + ur - ul)/(fld + frd);
            const auto cha = abs(ds - d0)/(0.5*(ds + d0));

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
            const auto c = sqrt(grav * ds);
            fl = 2.0 * (c-cl);
        }

        else
        {
            const auto ges = sqrt(0.5*grav*(ds + dl)/(ds*dl));
            fl = (ds - dl) * ges;
        }

        if (ds <= dr)
        {
            const auto c = sqrt(grav * ds);
            fr = 2.0 * (c-cr);
        }

        else
        {
            const auto ges = sqrt(0.5*grav*(ds + dr)/(ds*dr));
            fr = (ds - dr) * ges;
        }

        //exact Riemann solver sstar = ustar
        const auto us = 0.5*(ul + ur) + 0.5*(fr - fl); // u_star
        const auto cs = sqrt(grav*ds); // c_star

        /***********************  computation of u and d *******************/
        //left wave
        if (s <= us)
        {
            sol.velocityY = vl;
            //left shock
            if (ds >= dl)
            {
                const auto ql = sqrt((ds + dl)*ds/(2.0*dl*dl)); // flux left
                const auto sl = ul - cl*ql; // shock left

                //left side of shock
                if(s <= sl)
                {
                    sol.waterDepth = dl;
                    sol.velocityX = ul;
                }

                //right side of shock
                else
                {
                    sol.waterDepth = ds;
                    sol.velocityX = us;
                }
            }

            //left rarefaction
            else
            {
                const auto shl = ul-cl; //wave speed at head left

                //right side of rarefaction
                if (s <= shl)
                {
                    sol.waterDepth = dl;
                    sol.velocityX = ul;
                }

                else
                {
                    const auto stl = us -cs; //wave speed at tail left

                    //inside the rarefaction
                    if(s <= stl)
                    {
                        sol.velocityX = (ul + 2.0*cl + 2.0 * s)/3.0;
                        const auto c = (ul + 2.0* cl - s)/3.0;
                        sol.waterDepth = c*c/grav;
                    }

                    //inside the star region
                    else
                    {
                        sol.waterDepth = ds;
                        sol.velocityX = us;
                    }
                }
            }
        }

        //right wave
        else
        {
            sol.velocityY = vr;

            //right shock
            if(ds >= dr)
            {
                const auto qr = sqrt((ds + dr)*ds /(2.0*dr*dr)); // flux right
                const auto sr = ur + cr*qr; // shock right

                //right side of shock
                if(s >= sr)
                {
                    sol.waterDepth = dr;
                    sol.velocityX = ur;
                }

                //left side of shock
                else
                {
                    sol.waterDepth = ds;
                    sol.velocityX = us;
                }

            //right rarefaction
            }

            else
            {
                const auto shr  = ur + cr; //wave speed at head right

                //right side of Rarefaction
                if (s >= shr)
                {
                    sol.waterDepth = dr;
                    sol.velocityX = ur;
                }

                else
                {
                    const auto str = us + cs; //wave speed at tail right

                    //inside the rarefaction
                    if(s>=str)
                    {
                        sol.velocityX = (ur - 2.0*cr + 2.0*s)/3.0;
                        const auto c = (-ur + 2.0*cr + s)/3.0;
                        sol.waterDepth = c*c/grav;

                    }

                    //inside the star region
                    else
                    {
                        sol.waterDepth = ds;
                        sol.velocityX = us;
                    }
                }
            }
        }
    }
    //============================== Flux computation ===================

    sol.flux[0] = sol.waterDepth * sol.velocityX;
    sol.flux[1] = sol.waterDepth * sol.velocityX * sol.velocityX + 0.5 * grav * sol.waterDepth * sol.waterDepth;
    sol.flux[2] = sol.waterDepth * sol.velocityX * sol.velocityY;

    return sol;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
