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
#ifndef DUMUX_DISCRETIZATION_FICKIAN_DIFFUSION_COEFFICIENTS_HH
#define DUMUX_DISCRETIZATION_FICKIAN_DIFFUSION_COEFFICIENTS_HH

#include <cassert>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Container storing the diffusion coefficients required by Fick's law.
 *        Uses the minimal possible container size and provides unified access.
 * \tparam Scalar The type used for scalar values
 * \tparam numPhases Number of phases in the fluid composition
 * \tparam numComponents Number of components in the fluid composition
 * \tparam onlyTracers If false, this means that the main component of
 *                     a phase is part of the components. In this case,
 *                     the storage container is optimized with respect to
 *                     memory consumption as diffusion coefficients of the
 *                     main component of a phase in itself are not stored.
 *                     If true, all diffusion coefficients of all components
 *                     are stored
 */
template <class Scalar, int numPhases, int numComponents, bool onlyTracers = false>
class FickianDiffusionCoefficients;

//! General case (mpnc), for compositions containing the phases' main components
template <class Scalar, int numPhases, int numComponents>
class FickianDiffusionCoefficients<Scalar, numPhases, numComponents>
{
public:
    template<class DiffCoeffFunc>
    void update(DiffCoeffFunc& computeDiffCoeff)
    {
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
        {
            unsigned int compIIdx = phaseIdx;
            for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            {
                if (compIIdx == compJIdx)
                    continue;
                diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)] = computeDiffCoeff(phaseIdx, compIIdx, compJIdx);
            }
        }
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        std::cout << phaseIdx << " " << compIIdx << " " << compJIdx << " | " << getIndex_(phaseIdx, compIIdx, compJIdx) << " | ";
        sortComponentIndices_(compIIdx, compJIdx);
        std::cout << phaseIdx << " " << compIIdx << " " << compJIdx << " | " << getIndex_(phaseIdx, compIIdx, compJIdx) << " | ";

        return diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)];
    }


private:
    /*
     * \brief For the mpnc case, we have a rectangular diffCoeff matrix (numPhases x numComponents).
     *        As we do not track a diffusion coefficient for a component within it's phase,
     *        this numComponents is reduced by 1. (numPhases x (numComponents - 1))
    */
    std::array<Scalar, numPhases * (numComponents - 1)> diffCoeff_;

    /*
     * \brief This calculates the linear index of the diffusion coefficient matrix.
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

    void sortComponentIndices_(int& compIIdx, int& compJIdx) const
    { if (compIIdx > compJIdx) std::swap(compIIdx, compJIdx); }
};

//! Specialization for 1pnc & compositions containing the phases' main components
template <class Scalar, int numComponents>
class FickianDiffusionCoefficients<Scalar, 1, numComponents>
{
public:
    template<class DiffCoeffFunc>
    void update(DiffCoeffFunc& computeDiffCoeff)
    {
        for (unsigned int compJIdx = 1; compJIdx < numComponents; ++compJIdx)
        {
            diffCoeff_[getIndex_(0, 0, compJIdx)] = computeDiffCoeff(0, 0, compJIdx);
        }
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        assert(phaseIdx == 0);
        assert(phaseIdx != compJIdx);
        assert(compIIdx != compJIdx);
        assert(compIIdx < numComponents);
        assert(compJIdx < numComponents);
        assert(compIIdx <= compJIdx);
        sortComponentIndices_(compIIdx, compJIdx);
        return diffCoeff_[compJIdx-1];
    }

private:
    std::array<Scalar, numComponents-1> diffCoeff_;

    constexpr int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    { return  compJIdx-1; }

    void sortComponentIndices_(int& compIIdx, int& compJIdx) const
    { if (compIIdx > compJIdx) std::swap(compIIdx, compJIdx); }
};

//! Specialization for 2p2c & compositions containing the phases' main components
template <class Scalar>
class FickianDiffusionCoefficients<Scalar, 2, 2>
{
public:
    template<class DiffCoeffFunc>
    void update(DiffCoeffFunc& computeDiffCoeff)
    {
        diffCoeff_[getIndex_(0, 0, 1)] = computeDiffCoeff(0, 0, 1);
        diffCoeff_[getIndex_(1, 1, 0)] = computeDiffCoeff(1, 1, 0);
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        sortComponentIndices_(compIIdx, compJIdx);
        return diffCoeff_[phaseIdx];
    }

private:
    std::array<Scalar, 2> diffCoeff_;

    constexpr int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    { return  phaseIdx; }

    void sortComponentIndices_(int& compIIdx, int& compJIdx) const
    { if (compIIdx > compJIdx) std::swap(compIIdx, compJIdx); }
};

//! Specialization for 3p2c & compositions containing the phases' main components
template <class Scalar>
class FickianDiffusionCoefficients<Scalar, 3, 2>
{
public:
    template<class DiffCoeffFunc>
    void update(DiffCoeffFunc& computeDiffCoeff)
    {
        diffCoeff_[getIndex_(0, 0, 1)] = computeDiffCoeff(0, 0, 1);
        diffCoeff_[getIndex_(1, 1, 0)] = computeDiffCoeff(1, 1, 0);
        diffCoeff_[getIndex_(2, 1, 0)] = computeDiffCoeff(2, 1, 0);
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        sortComponentIndices_(compIIdx, compJIdx);
        return diffCoeff_[phaseIdx];
    }

private:
    std::array<Scalar, 3> diffCoeff_;

    constexpr int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    { return phaseIdx; }

    void sortComponentIndices_(int& compIIdx, int& compJIdx) const
    { if (compIIdx > compJIdx) std::swap(compIIdx, compJIdx); }
};

//! Specialization for mpnc & compositions that only contain tracers
template <class Scalar, int numPhases, int numComponents>
class FickianDiffusionCoefficients<Scalar, numPhases, numComponents, true>
{
public:
    template<class DiffCoeffFunc>
    void update(DiffCoeffFunc& computeDiffCoeff)
    {
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffCoeff_[getIndex_(phaseIdx, compIdx)]
                    = computeDiffCoeff(phaseIdx, phaseIdx, compIdx);
    }

    const Scalar& operator()(int phaseIdx, int compIIdx, int compJIdx) const
    {
        sortComponentIndices_(compIIdx, compJIdx);
        assert(phaseIdx < numPhases);
        assert(compIIdx < numComponents);
        assert(compJIdx < numComponents);
        assert(compIIdx <= compJIdx);
        return diffCoeff_[getIndex_(phaseIdx, compJIdx)];
    }

private:
    std::array<Scalar, numPhases*numComponents> diffCoeff_;

    constexpr int getIndex_(int phaseIdx, int compJIdx) const
    { return phaseIdx * numComponents + compJIdx; }

    void sortComponentIndices_(int& compIIdx, int& compJIdx) const
    { if (compIIdx > compJIdx) std::swap(compIIdx, compJIdx); }
};

} // end namespace Dumux

#endif
