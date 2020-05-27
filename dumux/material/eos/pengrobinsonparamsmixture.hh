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
 * \brief The mixing rule for the oil and the gas phases of the SPE5 problem.
 *
 * This problem comprises \f$H_2O\f$, \f$C_1\f$, \f$C_3\f$, \f$C_6\f$,
 * \f$C_10\f$, \f$C_15\f$ and \f$C_20\f$ as components.
 *
 * See:
 *
 * R. Reid, et al. (1987, pp. 43-44) \cite reid1987 <BR>
 *
 * and
 *
 * J.E. Killough, et al. (1987) \cite SPE5
 */
#ifndef DUMUX_PENG_ROBINSON_PARAMS_MIXTURE_HH
#define DUMUX_PENG_ROBINSON_PARAMS_MIXTURE_HH

#include <algorithm>
#include <dumux/material/constants.hh>

#include "pengrobinsonparams.hh"

namespace Dumux {

/*!
 * \ingroup EOS
 * \brief The mixing rule for the oil and the gas phases of the SPE5 problem.
 *
 * This problem comprises \f$H_2O\f$, \f$C_1\f$, \f$C_3\f$, \f$C_6\f$,
 * \f$C_10\f$, \f$C_15\f$ and \f$C_20\f$ as components.
 *
 * See:
 *
 * R. Reid, et al. (1987, pp. 43-44) \cite reid1987 <BR>
 *
 * and
 *
 * J.E. Killough, et al. (1987) \cite SPE5
 */
template <class Scalar, class FluidSystem, int phaseIdx, bool useSpe5Relations=false>
class PengRobinsonParamsMixture
    : public PengRobinsonParams<Scalar>
{
    enum { numComponents = FluidSystem::numComponents };

    // Peng-Robinson parameters for pure substances
    using PureParams = PengRobinsonParams<Scalar>;

public:
    /*!
     * \brief Update Peng-Robinson parameters for the pure components.
     * \param fluidState Thermodynamic state of the fluids
     *
     */
    template <class FluidState>
    void updatePure(const FluidState &fluidState)
    {
        updatePure(fluidState.temperature(phaseIdx),
                   fluidState.pressure(phaseIdx));
    }

    /*!
     * \brief Peng-Robinson parameters for the pure components.
     * \param temperature Temperature in \f$\mathrm{[K]}\f$
     * \param pressure pressure in \f$\mathrm{[Pa]}\f$
     * This method is given by the SPE5 paper.
     */
    void updatePure(Scalar temperature, Scalar pressure)
    {
        // Calculate the Peng-Robinson parameters of the pure
        // components
        //
        // See: R. Reid, et al.: The Properties of Gases and Liquids,
        // 4th edition, McGraw-Hill, 1987, p. 43
        for (int i = 0; i < numComponents; ++i) {
            Scalar pc = FluidSystem::criticalPressure(i);
            Scalar omega = FluidSystem::acentricFactor(i);
            Scalar Tr = temperature/FluidSystem::criticalTemperature(i);
            Scalar RTc = Constants<Scalar>::R*FluidSystem::criticalTemperature(i);

            Scalar f_omega;

            if (useSpe5Relations) {
                if (omega < 0.49) f_omega = 0.37464  + omega*(1.54226 + omega*(-0.26992));
                else              f_omega = 0.379642 + omega*(1.48503 + omega*(-0.164423 + omega*0.016666));
            }
            else
                f_omega = 0.37464 + omega*(1.54226 - omega*0.26992);

            using std::sqrt;
            Scalar tmp = 1 + f_omega*(1 - sqrt(Tr));
            tmp = tmp*tmp;

            Scalar a = 0.4572355*RTc*RTc/pc * tmp;
            Scalar b = 0.0777961 * RTc / pc;
            using std::isfinite;
            assert(isfinite(a));
            assert(isfinite(b));

            this->pureParams_[i].setA(a);
            this->pureParams_[i].setB(b);
        }

        updateACache_();
    }

    /*!
     * \brief Calculates the "a" and "b" Peng-Robinson parameters for
     *        the mixture.
     *
     * The updatePure() method needs to be called _before_ calling
     * this method!
     * \param fs the thermodynamic state of the fluids
     */
    template <class FluidState>
    void updateMix(const FluidState &fs)
    {
        Scalar sumx = 0.0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumx += fs.moleFraction(phaseIdx, compIdx);
        using std::max;
        sumx = max(1e-10, sumx);

        // Calculate the Peng-Robinson parameters of the mixture
        //
        // See: R. Reid, et al.: The Properties of Gases and Liquids,
        // 4th edition, McGraw-Hill, 1987, p. 82
        Scalar a = 0;
        Scalar b = 0;
        for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
            Scalar xi = fs.moleFraction(phaseIdx, compIIdx) / sumx;
            using::std::isfinite;
            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                Scalar xj = fs.moleFraction(phaseIdx, compJIdx) / sumx;

                // mixing rule from Reid, page 82
                a +=  xi * xj * aCache_[compIIdx][compJIdx];

                assert(isfinite(a));
            }

            // mixing rule from Reid, page 82
            b += xi * this->pureParams_[compIIdx].b();
            assert(isfinite(b));
        }

        this->setA(a);
        this->setB(b);
    }

    /*!
     * \brief Calculates the "a" and "b" Peng-Robinson parameters for
     *        the mixture provided that only a single mole fraction
     *        was changed.
     *
     * The updatePure() method needs to be called _before_ calling
     * this method!
     * \param fs the thermodynamic state of the fluids
     * \param compIdx the component index
     */
    template <class FluidState>
    void updateSingleMoleFraction(const FluidState &fs,
                                  int compIdx)
    {
        updateMix(fs);
    }

    /*!
     * \brief Return the Peng-Robinson parameters of a pure substance,
     * \param compIdx the component index
     */
    const PureParams &pureParams(int compIdx) const
    { return pureParams_[compIdx]; }

    /*!
     * \brief Returns the Peng-Robinson parameters for a pure component.
     * \param compIdx the component index
     */
    const PureParams &operator[](int compIdx) const
    {
        assert(0 <= compIdx && compIdx < numComponents);
        return pureParams_[compIdx];
    }

    [[deprecated("Will be removed after Dumux 3.3")]]
    void checkDefined() const {}

protected:
    PureParams pureParams_[numComponents];

private:
    void updateACache_()
    {
        for (int compIIdx = 0; compIIdx < numComponents; ++ compIIdx) {
            for (int compJIdx = 0; compJIdx < numComponents; ++ compJIdx) {
                // interaction coefficient as given in SPE5
                Scalar Psi = FluidSystem::interactionCoefficient(compIIdx, compJIdx);
                using std::sqrt;
                aCache_[compIIdx][compJIdx] =
                    sqrt(this->pureParams_[compIIdx].a()
                              * this->pureParams_[compJIdx].a())
                    * (1 - Psi);
            }
        }
    }

    Scalar aCache_[numComponents][numComponents];
};

} // end namespace

#endif
