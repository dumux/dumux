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
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for nitrogen and oxygen.
 */
#ifndef DUMUX_BINARY_COEFF_N2_O2_HH
#define DUMUX_BINARY_COEFF_N2_O2_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>

#include <dumux/material/components/o2.hh>
#include <dumux/material/components/n2.hh>

namespace Dumux
{
namespace BinaryCoeff
{

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for nitrogen and oxygen.
 */
class N2_O2
{
public:
    /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for molecular oxygen in liquid nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        DUNE_THROW(Dune::NotImplemented, "henry coefficient for oxygen in liquid nitrogen");
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular oxygen in liquid nitrogen.
     *
     * Uses Fuller's method to determine the diffusion of water in nitrogen.
     * This is only valid at "low" pressures.
     * See: R. Reid, et al. (1987, pp. 587-588) \cite reid1987
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        using N2 = Dumux::N2<Scalar>;
        using O2 = Dumux::O2<Scalar>;

        // atomic diffusion volumes
        static const std::array<Scalar, 2> SigmaNu = { 18.5 /* N2 */,  16.3 /* O2 */ };
        // molar masses [g/mol]
        static const std::array<Scalar, 2> M = { N2::molarMass()*1e3, O2::molarMass()*1e3 };

        static const Scalar Mab = harmonicMean(M[0], M[1]);

        using std::pow;
        using std::sqrt;
        using std::cbrt;
        static const Scalar tmp = cbrt(SigmaNu[0]) + cbrt(SigmaNu[1]);
        return 1e-4 * (143.0*pow(temperature, 1.75))/(pressure*sqrt(Mab)*tmp*tmp);
    }

    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular oxygen in liquid nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "diffusion coefficient for liquid oxygen and nitrogen");
    }
};

}
} // end namespace

#endif
