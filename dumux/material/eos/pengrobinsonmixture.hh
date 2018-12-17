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
 * \ingroup EOS
 * \brief Implements the Peng-Robinson equation of state for a
 *        mixture.
 */
#ifndef DUMUX_PENG_ROBINSON_MIXTURE_HH
#define DUMUX_PENG_ROBINSON_MIXTURE_HH

#include "pengrobinson.hh"

#include <dumux/material/constants.hh>

namespace Dumux {

/*!
 * \ingroup EOS
 * \brief Implements the Peng-Robinson equation of state for a
 *        mixture.
 */
template <class Scalar, class StaticParameters>
class PengRobinsonMixture
{
    enum { numComponents = StaticParameters::numComponents };
    using PengRobinson = Dumux::PengRobinson<Scalar>;

    // this class cannot be instantiated!
    PengRobinsonMixture() {};

    // the u and w parameters as given by the Peng-Robinson EOS
    static const Scalar u;
    static const Scalar w;

public:
    /*!
     * \brief Computes molar volumes \f$\mathrm{[m^3 / mol]}\f$ where the Peng-Robinson EOS is
     *        true.
     * \param Vm Molar Volume \f$\mathrm{[m^3 / mol]}\f$
     * \param fs Thermodynamic state of the fluids
     * \param params Parameters
     * \param phaseIdx The phase index
     * \return Number of solutions.
     */
    template <class MutableParams, class FluidState>
    static int computeMolarVolumes(Scalar *Vm,
                                   const MutableParams &params,
                                   int phaseIdx,
                                   const FluidState &fs)
    {
        return PengRobinson::computeMolarVolumes(Vm, params, phaseIdx, fs);
    }

    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in the phase.
     * \param fs Thermodynamic state of the fluids
     * \param params Parameters
     * \param phaseIdx The phase index
     * \param compIdx The index of the component
     *
     * The fugacity coefficient \f$\phi_i\f$ of a component \f$i\f$ is
     * defined as
     * \f[f_i = \phi_i x_i \;,\f]
     * where \f$f_i\f$ is the component's fugacity and \f$x_i\f$ is
     * the component's mole fraction.
     *
     * See:
     *
     * R. Reid, et al. (1987, pp. 42-44, 143-145) \cite reid1987
     */
    template <class FluidState, class Params>
    static Scalar computeFugacityCoefficient(const FluidState &fs,
                                             const Params &params,
                                             int phaseIdx,
                                             int compIdx)
    {
        // note that we normalize the component mole fractions, so
        // that their sum is 100%. This increases numerical stability
        // considerably if the fluid state is not physical.
        Scalar Vm = params.molarVolume(phaseIdx);

        // Calculate b_i / b
        Scalar bi_b = params.bPure(phaseIdx, compIdx) / params.b(phaseIdx);

        // Calculate the compressibility factor
        Scalar RT = Constants<Scalar>::R*fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx); // molar volume in [bar]
        Scalar Z = p*Vm/RT; // compressibility factor

        // Calculate A^* and B^* (see: Reid, p. 42)
        Scalar Astar = params.a(phaseIdx)*p/(RT*RT);
        Scalar Bstar = params.b(phaseIdx)*p/(RT);

        // calculate delta_i (see: Reid, p. 145)
        Scalar sumMoleFractions = 0.0;
        for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            sumMoleFractions += fs.moleFraction(phaseIdx, compJIdx);

        using std::sqrt;
        Scalar deltai = 2*sqrt(params.aPure(phaseIdx, compIdx))/params.a(phaseIdx);
        Scalar tmp = 0;
        for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
            tmp +=
                fs.moleFraction(phaseIdx, compJIdx)
                / sumMoleFractions
                * sqrt(params.aPure(phaseIdx, compJIdx))
                * (1.0 - StaticParameters::interactionCoefficient(compIdx, compJIdx));
        }
        deltai *= tmp;

        Scalar base =
            (2*Z + Bstar*(u + sqrt(u*u - 4*w))) /
            (2*Z + Bstar*(u - sqrt(u*u - 4*w)));
        Scalar expo =  Astar/(Bstar*sqrt(u*u - 4*w))*(bi_b - deltai);

        using std::exp;
        using std::max;
        using std::min;
        using std::pow;
        Scalar fugCoeff =
            exp(bi_b*(Z - 1))/max(1e-9, Z - Bstar) *
            pow(base, expo);

        ////////
        // limit the fugacity coefficient to a reasonable range:
        //
        // on one side, we want the mole fraction to be at
        // least 10^-3 if the fugacity is at the current pressure
        //

        fugCoeff = min(1e10, fugCoeff);
        //
        // on the other hand, if the mole fraction of the component is 100%, we want the
        // fugacity to be at least 10^-3 Pa
        //
        fugCoeff = max(1e-10, fugCoeff);
        ///////////
        using std::isfinite;
        if (!isfinite(fugCoeff)) {
            std::cout << "Non finite phi: " << fugCoeff << "\n";
        }

        return fugCoeff;
    }

};

template<class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::u = 2.0;
template<class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::w = -1.0;

} // end namespace Dumux

#endif
