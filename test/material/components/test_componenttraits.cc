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
 * \ingroup MaterialTests
 * \brief Test the component traits.
 */

#include "config.h"

#include <type_traits>

#include <dumux/material/components/air.hh>
#include <dumux/material/components/componenttraits.hh>

int main(int argc, char *argv[])
{
    using namespace Dumux;

    using Traits = ComponentTraits<Components::Air<double>>;
    static_assert(Traits::hasGasState, "Air component is reported to have no gas state?!");
    static_assert(!Traits::hasSolidState, "Air component is reported to implement a solid state?!");
    static_assert(!Traits::hasLiquidState, "Air component is reported to implement a liquid state?!");
    static_assert(!Traits::isIon, "Air component is reported to be an ion?!");
    static_assert(std::is_same<double, Traits::Scalar>::value, "Scalar type not correctly reported!");
}
