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
 * \brief Container storing the diffusion coefficients required by Fick's law.
 *        Uses the minimal possible container size and provides unified access.
 */
#ifndef DUMUX_FLUX_FICKIAN_DIFFUSION_COEFFICIENTS_HH
#define DUMUX_FLUX_FICKIAN_DIFFUSION_COEFFICIENTS_HH

#include <array>
#include <cassert>
#include <algorithm>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Container storing the diffusion coefficients required by Fick's law.
 *        Uses the minimal possible container size and provides unified access.
 * \tparam Scalar The type used for scalar values
 * \tparam numPhases Number of phases in the fluid composition
 * \tparam numComponents Number of components in the fluid composition
 */
template <class Scalar, int numPhases, int numComponents>
class FickianDiffusionCoefficients
{
public:
    template<class DiffCoeffFunc>
    void update(const DiffCoeffFunc& computeDiffCoeff)
    {
        // fill the diffusion coefficient, only compute the ones we need, see getIndex_ doc
        // the first component index "compIIdx" is always the main compnent index (phaseIdx)
        // if we have less components than phases we need to limit the index
        // this last case only occurs for 3p2c, there are no other special cases (2p1c doesn't have diffusion, max. number of phases is 3)
        static_assert(numPhases <= numComponents || (numPhases == 3 && numComponents == 2),
                      "This combination of numPhases and numComponents is not supported!");
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (int compIIdx = std::min(phaseIdx, numComponents-1), compJIdx = 0; compJIdx < numComponents; ++compJIdx)
                if (compIIdx != compJIdx)
                    diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)] = computeDiffCoeff(phaseIdx, compIIdx, compJIdx);
    }

    Scalar operator() (int phaseIdx, int compIIdx, int compJIdx) const
    {
        sortComponentIndices_(phaseIdx, compIIdx, compJIdx);
        assert(compIIdx != compJIdx);
        return diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)];
    }


private:
    /*
     * \brief For the mpnc case, we do not track a diffusion coefficient for the component of a phase.
     *        Hence, we have a matrix of dimensions (numPhases * (numComponents - 1))
     *        We store this matrix in a one-dimensional array
    */
    std::array<Scalar, numPhases * (numComponents - 1)> diffCoeff_;

    /*
     * \brief This calculates the array index to a corresponding matrix entry.
     *        An example with a 3p5c model:
     *                        compJIdx
     *                   [ x  0  1  2  3  ]
     *          phaseIdx [ 4  x  5  6  7  ]
     *                   [ 8  9  x  10 11 ]
     *        -- "(phaseIdx * (numComponents-1))" advances to the phaseIdx row (numComp-1 per phase)
     *        -- "+ compJIdx" advances to the component index
     *        -- "- static_cast<int>(phaseIdx < compJIdx)" if the index is after the diagonal, -1.
     */
    constexpr int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    { return (phaseIdx * (numComponents-1)) + compJIdx - static_cast<int>(phaseIdx < compJIdx); }

    void sortComponentIndices_(int phaseIdx, int& compIIdx, int& compJIdx) const
    { if (compIIdx != std::min(phaseIdx, numComponents-1)) std::swap(compIIdx, compJIdx); }
};

} // end namespace Dumux

#endif
