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
 * \brief Analytical solution of Mandel problem
 */
#ifndef DUMUX_MANDEL_ANALYTICAL_HH
#define DUMUX_MANDEL_ANALYTICAL_HH

#include <dune/common/fvector.hh>

namespace Dumux {

using Scalar = double;

template <class Scalar>
class MandelAnalyticalSolution{

public:
Scalar a ;
private:

/*!
 * \brief coefficient A1 for solution
 *
 * \return Scalar
 */
Scalar calculateA1()
{
    return 2.0 * biotCoefficient_ +
           2.0 * (lame_ + G_)/(biotModulus_ * biotCoefficient_);
}


Scalar calculateA2()
{
    return 2.0 * biotCoefficient_ * G_ / (lame_ + 2.0 * G);
}

Scalar force_; // force appied on the surface [N]
Scalar a_;      // width/2 [m]
Scalar b_;      // height/2 [m]

Dune::FieldVector<Scalar,3> biotCoefficient_;

Scalar biotModulus_; // [Pa]

Scalar mu_; // Poisson's ratio
Scalar E_;  // Young's modulus [Pa]

Scalar G_;
Scalar lame_;

Scalar Mxx;
Scalar Mxz;
Scalar Mzz;
};

}
#endif
