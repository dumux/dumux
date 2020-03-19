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
 * \ingroup FluxTests
 * \brief This test makes sure that the diffusion coefficient containers
 *        store the proper information, and are accessed correctly.
 */

#include <config.h>
#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>

// include all diffusioncoefficient containers
#include <dumux/flux/fickiandiffusioncoefficients.hh>
#include <dumux/flux/maxwellstefandiffusioncoefficients.hh>

template<int nP, int nC>
struct TestFluidSystem
{
     static constexpr int numPhases = nP;
     static constexpr int numComponents = nC;

     int binaryDiffusionCoefficient(const int phaseIdx, const int compIIdx, const int compJIdx) const
     { return phaseIdx*10000 + compIIdx*100 + compJIdx; }
};

int main()
{
    using namespace Dumux;
    using Scalar = int;

    const int numPhases = 3;
    const int numComponents = 8;

    TestFluidSystem<numPhases,numComponents> fluidSystem;

    // temporary for debugging
    // {
    for (int phaseIdx = 0; phaseIdx < fluidSystem.numPhases; ++ phaseIdx)
    {
        int compIIdx = phaseIdx;
        for (int compJIdx = 0; compJIdx < fluidSystem.numComponents; ++compJIdx)
        {
            std::cout << fluidSystem.binaryDiffusionCoefficient(phaseIdx,compIIdx,compJIdx) << "  ";
        }
        std::cout << "\n";
    }
    // }

    const bool onlyTracers = false;
    using DiffusionCoefficientsContainer = FickianDiffusionCoefficients<Scalar,
                                                                        numPhases,
                                                                        numComponents,
                                                                        onlyTracers>;

    auto getDiffusionCoefficient = [&](int phaseIdx, int compIIdx, int compJIdx)
    {
        return fluidSystem.binaryDiffusionCoefficient(phaseIdx,
                                                      compIIdx,
                                                      compJIdx);
    };

    DiffusionCoefficientsContainer mpncDiffCoeffs;
    mpncDiffCoeffs.update(getDiffusionCoefficient);

    for (int phaseIdx = 0; phaseIdx < fluidSystem.numPhases; ++ phaseIdx)
    {
        int compIIdx = phaseIdx;
        std::cout << phaseIdx << " : \n";
        for (int compJIdx = 0; compJIdx < fluidSystem.numComponents; ++compJIdx)
        {
            if (compIIdx == compJIdx)
                continue;
            int binaryCoeff = fluidSystem.binaryDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
            int containerCoeff = mpncDiffCoeffs(phaseIdx, compIIdx, compJIdx);
            if (containerCoeff != binaryCoeff)
                DUNE_THROW(Dune::InvalidStateException, "Coefficients differ! " << "fluidSystem: " << binaryCoeff << " and "
                                                                                << "container: " << containerCoeff <<
                                                        ", for indicies (p,cI,cJ): (" << phaseIdx << ","
                                                                                      << compIIdx << ","
                                                                                      << compJIdx<<")" );
            std::cout << " \n";
        }
    }

    for (int phaseIdx = 0; phaseIdx < fluidSystem.numPhases; ++ phaseIdx)
    {
        int compJIdx = phaseIdx;
        std::cout << phaseIdx << " : \n";
        for (int compIIdx = 0; compIIdx < fluidSystem.numComponents; ++compIIdx)
        {
            if (compJIdx == compIIdx)
                continue;
            int binaryCoeff = fluidSystem.binaryDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
            int containerCoeff = mpncDiffCoeffs(phaseIdx, compIIdx, compJIdx);
            if (containerCoeff != binaryCoeff)
                DUNE_THROW(Dune::InvalidStateException, "Coefficients differ! " << "fluidSystem: " << binaryCoeff << " and "
                                                                                << "container: " << containerCoeff <<
                                                        ", for indicies (p,cI,cJ): (" << phaseIdx << ","
                                                                                      << compIIdx << ","
                                                                                      << compJIdx<<")" );
            std::cout << " \n";
        }
    }


//     // test 1pnc
//     template <class Scalar, int numComponents>
//     class FickianDiffusionCoefficients<Scalar, 1, numComponents>

//     // test 2p2c
//     template <class Scalar>
//     class FickianDiffusionCoefficients<Scalar, 2, 2>

//     // test 3p2c
//     template <class Scalar>
//     class FickianDiffusionCoefficients<Scalar, 3, 2>

//     // test tracer mpnc
//     template <class Scalar, int numPhases, int numComponents>
//     class FickianDiffusionCoefficients<Scalar, numPhases, numComponents, true>

//     // test max stef
//     template <class Scalar, int numPhases, int numComponents>
//     class MaxwellStefanDiffusionCoefficients

    return 0.0;
}
