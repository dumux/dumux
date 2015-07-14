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
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_N2CONSTVISCOSITY_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_N2CONSTVISCOSITY_HH

/*!
 * \file
 * \ingroup Components
 *
 * \brief Properties of pure molecular nitrogen \f$N_2\f$ with constant viscosity.
 *
 * The constant viscosity is useful to get a desired Reynolds number.
 */

#include <dumux/material/components/n2.hh>

namespace Dumux
{
  /**
   * \ingroup BoxStokesModel
   * \ingroup ImplicitTestProblems
   *
   * \brief Properties of pure molecular nitrogen \f$N_2\f$ with constant
   *        viscosity.
   *
   * The constant viscosity is useful to get a desired Reynolds number.
   *
   * \tparam Scalar The type used for scalar values
   */
  template <class Scalar>
  class N2ConstViscosity : public N2<Scalar>
  {
    typedef Dumux::IdealGas<Scalar> IdealGas;

  public:
    /**
     * \brief A human readable name for nitrogen with fixed viscosity.
     */
    static const char *name()
    {
      return "N2constViscosity";
    }

    /**
     * \brief The dynamic viscosity given by a fixed value to guarantee
     *        a Reynolds number.
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
      Scalar reynoldsNumber = 100.0;
      Scalar characteristicLength = 1.0;
      Scalar characteristicVelocity = 1.0;
      // Density of N2 for pressure of 0.1 MPa,
      // temperature of 283.15 degree Celsius (from NIST homepage)
      Scalar characteristicDensity = 1.1903;

      return characteristicLength * characteristicVelocity * characteristicDensity
             / reynoldsNumber;
    }
  };
} // end namespace
#endif // DUMUX_TEST_FREEFLOW_NAVIERSTOKES_N2CONSTVISCOSITY_HH
