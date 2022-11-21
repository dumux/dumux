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
 * \ingroup PoromechanicsTests
 * \brief Index to detect the region
 */

#ifndef DUMUX_APP_2PINJECTION_REGION_INDEX_HH
#define DUMUX_APP_2PINJECTION_REGION_INDEX_HH

#include <cmath>
#include <map>
namespace Dumux{

enum class ZoneType{
    LeftFault,
    RightFault,
    Caprock,
    Reservoir
};

constexpr int ZoneNum = 4;
std::array<ZoneType,ZoneNum> AllZoneTypes = {ZoneType::LeftFault,ZoneType::RightFault,ZoneType::Caprock,ZoneType::Reservoir};

const std::map<ZoneType,std::string> ZoneName =
    {{ZoneType::LeftFault,"LeftFault"},
     {ZoneType::RightFault,"RightFault"},
     {ZoneType::Caprock,"Caprock"},
     {ZoneType::Reservoir,"Reservoir"}};

template<class GlobalPosition>
ZoneType zoneTypeAtPos(const GlobalPosition& globalPos)
{
    using std::tan;
    using std::abs;

    // check on Fault
    // Fault 1000m long
    // (x - (-1500))
    // each fault zone has width of 10 m
    const double width = 10.0;
    const GlobalPosition leftFaultCenter = {500,-1500};
    const GlobalPosition rightFaultCenter = {1500,-1500};

    if (abs(globalPos[1] - leftFaultCenter[1]) < 500)
    {
        static const double tan_80 = tan(80.0/180*M_PI);
        //std::cout << M_PI << std::endl;
        //std::cout << tan_80 << std::endl;
        double localLeftFractureCenter = (globalPos[1] - leftFaultCenter[1]) /tan_80 + leftFaultCenter[0];

        // the left fault
        if (std::abs(localLeftFractureCenter - globalPos[0]) < width/2.0)
            return ZoneType::LeftFault;

        // the right fault
        if (std::abs(localLeftFractureCenter + rightFaultCenter[0] - leftFaultCenter[0] - globalPos[0]) < width/2.0)
            return ZoneType::RightFault;
    }

    // the upper caprock
    if (globalPos[1] > -1450 && globalPos[1] < -1300)
        return ZoneType::Caprock;

    // the lower caprock
    if (globalPos[1] > -1700 && globalPos[1] < -1550)
        return ZoneType::Caprock;

    // else
    return ZoneType::Reservoir;
}


}
#endif
