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

#include <config.h>
#include <dune/common/exceptions.hh>
#include <dumux/geomechanics/stressstate/math/mohrspace.hh>

int main(int argc, char *argv[])
{
    using namespace Dumux;
    Point p(1.0,3.0);
    Line l(4.27,11.5);

    [[maybe_unused]] double dist = pointLineDistance<double>(p,l);
    [[maybe_unused]] double solution = 2.9118;
    if (std::abs(dist - solution) > 1e-3)
    {DUNE_THROW(Dune::MathError, "test failed.");}
}
