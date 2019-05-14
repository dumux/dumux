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
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Manning.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */
#ifndef DUMUX_FRICTIONLAW_MANNING_HH
#define DUMUX_FRICTIONLAW_MANNING_HH

#include <algorithm>
#include <cmath>
#include "nikuradse.hh"

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Manning.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename Scalar>
class FrictionLawManning : FrictionLawNikuradse<Scalar>
{
public:
    /*!
     * \brief Compute the friction ustar_h.
     *
     * \param h water depth.
     * \param manningN Mannings friction value.
     * \return ustar_h friction used for the source term in shallow water models.
     */

    Scalar computeUstarH(const Scalar h,const Scalar manningN, const Scalar gravity)
    {
        using std::pow;

        Scalar ustar_h = 0.0;
        Scalar rough_h = pow(25.68/(1.0/manningN),6.0);

        rough_h = this->limitRoughH(rough_h, h);

        auto cfric = pow((h + rough_h),1.0/6.0) * 1.0/(manningN);
        ustar_h = gravity / pow(cfric,2.0);
        return ustar_h;
    }

};

} // end namespace Dumux

#endif // DUMUX_FRICTIONLAW_MANNING_HH
