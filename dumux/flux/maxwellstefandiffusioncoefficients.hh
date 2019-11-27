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
#ifndef DUMUX_DISCRETIZATION_MAXWELLSTEFAN_DIFFUSION_COEFFICIENTS_HH
#define DUMUX_DISCRETIZATION_MAXWELLSTEFAN_DIFFUSION_COEFFICIENTS_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

// General case
template <class Scalar, int numPhases, int numComponents>
class MaxwellStefanDiffusionCoefficients
{
public:
    template<class DiffCoeffFunc>
    void update(DiffCoeffFunc& computeDiffCoeff)
    {
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (unsigned int compIIdx = 0; compIIdx < numComponents; ++compIIdx)
                for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
                    if(compIIdx != compJIdx && compIIdx < compJIdx)
                        diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)] = computeDiffCoeff(phaseIdx, compIIdx, compJIdx);
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        assert(phaseIdx != compJIdx);
        return diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)];
    }

private:
    std::array<Scalar, (numPhases * ((numComponents * numComponents - numComponents)/2))> diffCoeff_;
    const int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    {
        return  phaseIdx * ((numComponents * numComponents - numComponents) / 2)
              + compIIdx * numComponents
              - ((compIIdx * compIIdx + compIIdx) / 2)
              + compJIdx - (compIIdx +1);
    }
};

} // end namespace Dumux

#endif
