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
#include <dune/common/float_cmp.hh>
#include <dumux/geomechanics/stressstate/stressdroplaw.hh>
#include <dumux/geomechanics/stressstate/stressdroplawparams.hh>
int main(int argc, char *argv[])
{
    using namespace Dumux;

    using Tensor = Dune::FieldMatrix<double,2,2>;
    using MohrSpaceTypeTraits = MohrSpaceTypeTraits<double,Tensor>;
    using Point = typename MohrSpaceTypeTraits::PointType;
    using Line = typename MohrSpaceTypeTraits::LineType;
    using MohrCircle = typename MohrSpaceTypeTraits::MohrCircleType;

    using StressDropLaw = StressDropLaw<double, Tensor, MohrSpaceTypeTraits>;
    using StressDropLawParams = StressDropLawParams<double, Line>;
    Tensor stress;
    stress[0][0] = -1.0;
    stress[0][1] = 2.0;
    stress[1][0] = 2.0;
    stress[1][1] = -4.0;

    StressDropLaw law(stress);
    StressDropLaw law2(stress);
    StressDropLawParams param(/*angle*/-45.0,/*cohesion*/0.0, /*stressdrop*/1.0);
    StressDropLawParams param2(/*angle*/-45.0,/*cohesion*/1.02, /*stressdrop*/1.0);

    if (!law.hasFailure(param))
    {
        DUNE_THROW(Dune::Exception,"Failure check was wrong.");
    }

    using Dune::FloatCmp::ne;
    using Dune::FloatCmp::eq;

    // test stress drop raw w/o rotation
    law.doStressDrop(param);
    if(ne(law.stressTensor()[0][0],-1.0) || ne(law.stressTensor()[0][1],1.0))
    {
        DUNE_THROW(Dune::Exception,"Stress drop without rotation was wrong.");
    }

    // test calculation of angle
    Point p1(0.0,0.0), p2(-1.5,-1), p0(-2.5,0.0);
    if(ne(law.angleBtwVectors(p1,p2,p0),45.0/180*M_PI))
    {
        DUNE_THROW(Dune::Exception,"Rotation angle was wrong.");
    }

    // test stress drop law in lower half of the mohr circle
    law2.doStressDrop(param2);
    auto stress2 = law2.stressTensor();
    Tensor origStress;
    origStress[0][0] = -0.95187992;
    origStress[0][1] = 1.00115844;
    origStress[1][0] = 1.00115844;
    origStress[1][1] = -4.04812008;
    for(int i=0; i<2; ++i)
    {
        for(int j=0; j<2; ++j)
        {
            if(ne(stress2[i][j],origStress[i][j],1e-4))
            {
                DUNE_THROW(Dune::Exception,"Stress drop with rotation was wrong.");
            }
        }
    }

    // test stress drop raw in upper half of mohr circle
    Tensor stress3;
    stress3[0][0] = -1.0;
    stress3[0][1] = -2.0;
    stress3[1][0] = -2.0;
    stress3[1][1] = -5.0;

    StressDropLawParams param3(/*angle*/-65.0,/*cohesion*/0.2, /*stressdrop*/1.0);
    StressDropLaw law3(stress3);
    law3.doStressDrop(param3);
    Tensor origStress3;
    origStress3[0][0] = -0.78560085;
    origStress3[0][1] = -1.02325387;
    origStress3[1][0] = -1.02325387;
    origStress3[1][1] = -5.21439915;
    //std::cout << "Stress 3" << "\n"
    //          << law3.stressTensor() <<std::endl;
    for(int i=0; i<2; ++i)
    {
        for(int j=0; j<2; ++j)
        {
            if(ne(law3.stressTensor()[i][j],origStress3[i][j],1e-3))
            {
                DUNE_THROW(Dune::Exception,"Stress drop with rotation was wrong.");
            }
        }
    }
    return 0;
}