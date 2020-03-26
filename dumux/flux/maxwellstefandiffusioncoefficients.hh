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
 * \ingroup Flux
 * \brief Container storing the diffusion coefficients required by the Maxwell-
 *        Stefan diffusion law. Uses the minimal possible container size and
 *        provides unified access.
 */
#ifndef DUMUX_FLUX_MAXWELLSTEFAN_DIFFUSION_COEFFICIENTS_HH
#define DUMUX_FLUX_MAXWELLSTEFAN_DIFFUSION_COEFFICIENTS_HH

#include <array>
#include <cassert>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Container storing the diffusion coefficients required by the Maxwell-
 *        Stefan diffusion law. Uses the minimal possible container size and
 *        provides unified access.
 * \tparam Scalar The type used for scalar values
 * \tparam numPhases Number of phases in the fluid composition
 * \tparam numComponents Number of components in the fluid composition
 */
template <class Scalar, int numPhases, int numComponents>
class MaxwellStefanDiffusionCoefficients
{
public:
    template<class DiffCoeffFunc>
    void update(const DiffCoeffFunc& computeDiffCoeff)
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx)
                for (int compJIdx = compIIdx+1; compJIdx < numComponents; ++compJIdx)
                    diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)]
                        = computeDiffCoeff(phaseIdx, compIIdx, compJIdx);
    }

    Scalar operator() (int phaseIdx, int compIIdx, int compJIdx) const
    {
        sortComponentIndices_(compIIdx, compJIdx);
        assert(compIIdx != compJIdx);
        return diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)];
    }

private:
    /*!
     * \brief Maxwell Stefan diffusion coefficient container.
     *        This container is sized to hold all required diffusion coefficients.
     *
     *        For each phase "(numPhases * ("
     *        We have a square coefficient matrix sized by
     *        the number of components "((numComponents * numComponents)".
     *        The diagonal is not used and removed " - numComponents)".
     *        The matrix is symmetrical, but only the upper triangle is required " / 2))".
     */
    std::array<Scalar, (numPhases * ((numComponents * (numComponents - 1)) / 2))> diffCoeff_;

    /*!
     * \brief Index logic for collecting the correct diffusion coefficient from the container.
     *
     *        First, we advance our index to the correct phase coefficient matrix:
     *        " phaseIdx * ((numComponents * numComponents - numComponents) / 2) ".
     *        The individual index within each phase matrix is then calculated using
     *        " i*n - (i^2+i)/2 + j-(i+1) ".
     *
     *        This index calculation can be reduced from the following:
     *        https://stackoverflow.com/questions/27086195/
     */
    constexpr int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    {
        return phaseIdx * ((numComponents * (numComponents - 1)) / 2)
               + compIIdx * numComponents
               - ((compIIdx * (compIIdx + 1)) / 2)
               + compJIdx - (compIIdx +1);
    }

    void sortComponentIndices_(int& compIIdx, int& compJIdx) const
    { if (compIIdx > compJIdx) std::swap(compIIdx, compJIdx); }
};

} // end namespace Dumux

#endif
