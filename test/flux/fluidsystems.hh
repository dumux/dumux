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
 * \brief Fluid systems for testing the diffusion containers
 */
#ifndef DUMUX_FLUIDSYSTEMS_HH
#define DUMUX_FLUIDSYSTEMS_HH

#include <string>
#include <vector>

namespace Dumux {

/*!
 * \brief A Fluid System for testing the MpNc container
 */
class MpNcFluidSystem
{
public:
    //! Constructor
    MpNcFluidSystem(const int numPhases, const int numComponents)
    : numPhases_(numPhases)
    , numComponents_(numComponents)
    {
        fillDiffCoeffMatrix_();
    }

    int binaryDiffusionCoefficient(const int phaseIdx, const int compIIdx, const int compJIdx) const
    { return binaryDiffusionCoefficients_[phaseIdx][compJIdx]; }

    int numPhases() const
    { return numPhases_; }

    int numComponents() const
    { return numComponents_; }

private:
    void fillDiffCoeffMatrix_()
    {
        int count = 0;
        binaryDiffusionCoefficients_.resize(numPhases_, std::vector<int>(numComponents_, -1));
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases_; ++ phaseIdx)
        {
            for (unsigned int compJIdx = 0; compJIdx < numComponents_; ++compJIdx)
            {
                if (compJIdx == phaseIdx)
                    continue;
                else
                {
                    binaryDiffusionCoefficients_[phaseIdx][compJIdx] = count;
                    count++;
                }
            }
        }
    }

    int numPhases_;
    int numComponents_;
    std::vector<std::vector<int>> binaryDiffusionCoefficients_;
};

} // end namespace Dumux

#endif
