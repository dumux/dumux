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
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dumux/geomechanics/stressstate/math/mohrspace.hh>

int main(int argc, char *argv[])
{
    using namespace Dumux;
    using Tensor = Dune::FieldMatrix<double,2,2>;
    using Point = Point<double>;
    using MohrCircle = MohrCircle<double, Point, Tensor>;
    using Line = Line<double>;

    Tensor stress;
    stress[0][0] = -1.0;
    stress[0][1] = 2.0;
    stress[1][0] = 2.0;
    stress[1][1] = -4.0;

    MohrCircle mohrCircle(stress);

    Line l(-0.3,0.7);

    using std::abs;
    if (abs(mohrCircle.center().x() + 2.5) > 1e-4){
        DUNE_THROW(Dune::Exception,"circle center is wrong.");
    }

    if (abs(mohrCircle.radius() - 2.5) > 1e-4){
        DUNE_THROW(Dune::Exception,"circle radius is wrong.");
    }

    if (!mohrCircle.hasIntersection(l)){
        DUNE_THROW(Dune::Exception,"intersection check is wrong.");
    }

    auto intersect = mohrCircle.intersections(l);
    std::array<Point,2> solution;
    solution[0] = Point(-0.1099,0.733);
    solution[1] = Point(-4.092,1.9276);
    for(int i =0; i<2 ; ++i)
    {
        //std::cout<< intersect[i].x() << " " << intersect[i].y() << std::endl;
        // dist has the type as dune::fieldvector<2>
        auto dist = intersect[i].pos() - solution[i].pos();
        if (dist.two_norm() > 1e-3){
        DUNE_THROW(Dune::Exception,"intersection is wrong.");
        }
    }
}
