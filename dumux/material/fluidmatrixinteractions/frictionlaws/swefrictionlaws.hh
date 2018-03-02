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
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of different friction laws for the shallow water
 *        equations. This class bundles all friction laws into one file, in
 *        practical applications, the friction value(s) are stored at the
 *        nodes/cell centers. One can also use different friction laws in
 *        one application.
 *
 */
#ifndef DUMUX_SWE_FRICTIONLAWS_HH
#define DUMUX_SWE_FRICTIONLAWS_HH
#include <algorithm>
#include <cmath>

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 *
 * \brief Implementation of different friction laws for the shallow water
 *        equations.
 *
 *
 */
template <class TypeTag>
class SweFrictionlaws
{

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    static const int manning = GET_PROP_VALUE(TypeTag,Manning);
    static const int chezy = GET_PROP_VALUE(TypeTag,Chezy);
    static const int nikuradse = GET_PROP_VALUE(TypeTag,Nikuradse);

public:


    /*!
     * \brief The computation of ustar_h for different friction laws.
     *
     * TODO describe all friction laws formulas
     *
     * \param h     water depth
     * \param ks    friction value, depending on the law
     * \params grav gravitation
     * \param law   integer id, defining the friction law to use
     * \return ustar_h.
     *
     */
    static Scalar ustar_h(const Scalar ks, const Scalar h, const Scalar grav, const int law)
    {
        Scalar ustar_h;
        using std::pow;
        using std::max;
        using std::min;
        using std::pow;
        using std::log;


        //compute friction
        if (law == manning){
            //Manning friction
            ustar_h = manningUstarH(ks,h,grav);
        }else if(law == chezy){
            //Chezy friction law
            ustar_h = grav / pow(ks,2.0);
        }else if(law == nikuradse){
            //Nikuradse friction
            //Nikuradse law uses a rough_h approximation for small depths to avoid extreme frictions
            Scalar mobility = 1.0;
            Scalar krw = 1.0;
            Scalar sw = 0.0;
            Scalar minUpperH;
            Scalar rough_h;
            const Scalar letL = 0.0, letT = 2.0, letE = 1.0;
            rough_h = ks;
            minUpperH = rough_h* 2.0; //we start to limit the friction for small water depths

            //Mobility after LET-model
            sw = min(h * (1.0/minUpperH),1.0);
            sw = max(0.0,sw);
            mobility = (krw * pow(sw,letL))/(pow(sw,letL) + letE * pow(1.0-sw,letT));
            rough_h = rough_h * (1.0 - mobility);

            //Nikuradse law
            ustar_h = pow(0.41,2.0)/pow(log((12*(h+ rough_h))/ks),2.0);
        }else{
            //now law is used and ustar_h is zero
            ustar_h = 0.0;
        }

        return ustar_h;
    }
};

} // end namespace Dumux

#endif // SWE_FRICTIONLAWS_HH
