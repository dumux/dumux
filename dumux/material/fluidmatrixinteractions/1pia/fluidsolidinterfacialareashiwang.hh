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
 * \brief Description of a interfacial area between solid and fluid phase
 */
#ifndef FLUIDSOLID_INTERFACIALAREA_SHI_WANG_HH
#define FLUIDSOLID_INTERFACIALAREA_SHI_WANG_HH

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Description of a interfacial area between solid and fluid phase
 */
template<class Scalar>
class FluidSolidInterfacialAreaShiWang
{
public:
    /*!
     * \brief Relation for the interfacial area between a fluid and a solid phase
     * after Shi & Wang, Transport in porous media (2011)
     *
     * \return interfacial area
     */
    static Scalar fluidSolidInterfacialArea(const Scalar porosity,
                                            const Scalar characteristicLength)
    { return 6.0*(1.0-porosity)/characteristicLength; }
};

} // end namespace Dumux

#endif
