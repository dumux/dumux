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
 * \ingroup OnePNCMinTests
 * \brief Corrected material properties of pure Calcium-Oxide \f$CaO\f$ without considering a porosity
 * change in the reaction of Calciumoxyde and Calciumhydroxyde.
 */

#ifndef DUMUX_MODIFIED_CAO_HH
#define DUMUX_MODIFIED_CAO_HH


#include <dumux/material/components/cao.hh>

namespace Dumux {
namespace Components {
/*!
 * \ingroup OnePNCMinTests
 * \brief A class for the ModifiedCaO properties.
 *
 * This class uses a different CaO density. It is to be used for calculating the
 * chemical reaction of CaO to Ca(OH)2 without considering the porosity change
 * according to Shao et al. (2013) \cite shao2013.
 */
template <class Scalar>
class ModifiedCaO : public  Components::CaO<Scalar>
{
public:

    /*!
     * \brief The corrected mass density \f$\mathrm{[kg/m^3]}\f$ of CaO.
     *
     * This density is to be used for calculating the chemical reaction
     * of CaO to Ca(OH)2 without considering the solid volume change.
     * See Shao et al. (2013) \cite shao2013.
     */
    static Scalar solidDensity(Scalar temperature)
    {
        return 1656;
    }

};
} // end namespace Components
} // end namespace Dumux

#endif
