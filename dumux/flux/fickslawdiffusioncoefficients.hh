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
#ifndef DUMUX_DISCRETIZATION_FICKS_LAW_DIFFUSION_COEFFICIENTS_HH
#define DUMUX_DISCRETIZATION_FICKS_LAW_DIFFUSION_COEFFICIENTS_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

// forward declaration
template <class Scalar, int numPhases, int numComponents>
class FicksLawDiffusionCoefficients
{
public:
    Scalar diffusionCoefficient(int phaseIdx, int compIdxI, int compIdxJ) const
    {
        if (compIdxJ < phaseIdx)
            return diffCoeff_[phaseIdx][compIdxJ];
        else if (compIdxJ > phaseIdx)
            return diffCoeff_[phaseIdx][compIdxJ-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient called for phaseIdx = compIdxJ");
    }

private:
    std::array<std::array<Scalar, numComponents-1>, 2> diffCoeff_;
};

template <class Scalar, int numComponents>
class FicksLawDiffusionCoefficients<Scalar, 1, numComponents>
{
public:
    Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx)
    {
        assert(phaseIdx == 0);
        assert(phaseIdx != compJIdx);
        assert(compIIdx != compJIdx);
        assert(compIIdx < numComponents);
        assert(compJIdx < numComponents);
        assert(compIIdx < compJIdx);
        return diffCoeff_[compJIdx-1];
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        assert(phaseIdx == 0);
        assert(phaseIdx != compJIdx);
        assert(compIIdx != compJIdx);
        assert(compIIdx < numComponents);
        assert(compJIdx < numComponents);
        assert(compIIdx < compJIdx);
        return diffCoeff_[compJIdx-1];
    }

private:
    std::array<Scalar, numComponents-1> diffCoeff_;
};

template <class Scalar>
class FicksLawDiffusionCoefficients<Scalar, 2, 2>
{
public:

    Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx)
    {
        assert(phaseIdx != compJIdx);
        return diffCoeff_[phaseIdx];
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        assert(phaseIdx != compJIdx);
        return diffCoeff_[phaseIdx];
    }

private:
     std::array<Scalar, 2> diffCoeff_;
};




} // end namespace Dumux

#endif
