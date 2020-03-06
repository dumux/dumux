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
 * \tparam allTracerComponents If false, this means that the main component of
 *                             a phase is part of the components. In this case,
 *                             the storage container is optimized with respect to
 *                             memory consumption as diffusion coefficients of the
 *                             main component of a phase in itself are not stored.
 *                             If true, all diffusion coefficients of all components
 *                             are stored
 */
template <class Scalar, int numPhases, int numComponents, bool allTracerComponents>
class FickianDiffusionCoefficients;

//! General case (mpnc), for compositions containing the phases' main components
template <class Scalar, int numPhases, int numComponents>
class FickianDiffusionCoefficients<Scalar, numPhases, numComponents, false>
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
        assert(phaseIdx != compJIdx);
        return diffCoeff_[getIndex_(phaseIdx, compIIdx, compJIdx)];
    }

private:
    std::array<Scalar, numPhases * (numComponents - 1)> diffCoeff_;
    int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    {
        return ((phaseIdx * (numComponents-1)) + compJIdx) - (phaseIdx < compJIdx);
    }
};

//! Specialization for 1pnc & compositions containing the phases' main components
template <class Scalar, int numComponents>
class FickianDiffusionCoefficients<Scalar, 1, numComponents, false>
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
        assert(compIIdx < compJIdx);
        return diffCoeff_[compJIdx-1];
    }

private:
    std::array<Scalar, numComponents-1> diffCoeff_;
    int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    { return  compJIdx-1; }
};

//! Specialization for 2p2c & compositions containing the phases' main components
template <class Scalar>
class FickianDiffusionCoefficients<Scalar, 2, 2, false>
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
        assert(phaseIdx != compJIdx);
        return diffCoeff_[phaseIdx];
    }

private:
    std::array<Scalar, 2> diffCoeff_;
    int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    { return  phaseIdx; }
};

//! Specialization for 3p2c & compositions containing the phases' main components
template <class Scalar>
class FickianDiffusionCoefficients<Scalar, 3, 2, false>
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
        assert(phaseIdx != compJIdx);
        return diffCoeff_[phaseIdx];
    }

private:
    std::array<Scalar, 3> diffCoeff_;
    int getIndex_(int phaseIdx, int compIIdx, int compJIdx) const
    { return  phaseIdx; }
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
        assert(phaseIdx < numPhases);
        assert(compIIdx != compJIdx);
        assert(compIIdx < numComponents);
        assert(compJIdx < numComponents);
        assert(compIIdx < compJIdx);
        return diffCoeff_[getIndex_(phaseIdx, compJIdx)];
    }

private:
    std::array<Scalar, numPhases*numComponents> diffCoeff_;
    int getIndex_(int phaseIdx, int compJIdx) const
    { return phaseIdx*numComponents + compJIdx; }
};

} // end namespace Dumux

#endif
