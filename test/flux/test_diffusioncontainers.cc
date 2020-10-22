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
#include <iostream>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>

// include all diffusioncoefficient containers
#include <dumux/flux/fickiandiffusioncoefficients.hh>
#include <dumux/flux/maxwellstefandiffusioncoefficients.hh>

template<int nP, int nC>
struct TestFluidSystem
{
    static constexpr int numPhases = nP;
    static constexpr int numComponents = nC;

    int binaryDiffusionCoefficient(int phaseIdx, int compIIdx, int compJIdx) const
    {
        if (compIIdx < compJIdx) std::swap(compIIdx, compJIdx);
        return phaseIdx*10000 + compIIdx*100 + compJIdx;
    }
};

template<class Container, int numPhases, int numComponents>
void testFickianContainer()
{
    TestFluidSystem<numPhases, numComponents> fs;
    Container diffCoeff;

    std::cout << "Testing " << Dune::className(diffCoeff) << std::endl;

    // fill the container with coefficients
    diffCoeff.update([fs](int phaseIdx, int compIIdx, int compJIdx) {
        return fs.binaryDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
    });

    // we test that the internal mapping is unique and that the interface is symmetric
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        for (int compIIdx = std::min(phaseIdx, numComponents-1), compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        {
            if (compIIdx != compJIdx)
            {
                const auto refD = fs.binaryDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
                const auto D = diffCoeff(phaseIdx, compIIdx, compJIdx);
                const auto D2 = diffCoeff(phaseIdx, compJIdx, compIIdx);
                if (D != refD)
                    DUNE_THROW(Dune::Exception,
                        "Retrieved incorrect diffusion coefficient for (p,c0,c1) = ("
                        << phaseIdx << "," << compIIdx << "," << compJIdx << "). "
                        << "Result: " << D << ", expected result: " << refD);
                if (D2 != D)
                    DUNE_THROW(Dune::Exception,
                        "Diffusion coefficient interface not symmetric in component indices for (p,c0,c1) = ("
                        << phaseIdx << "," << compIIdx << "," << compJIdx << "). "
                        << "Result: " << D2 << ", expected result: " << D);
            }
        }
    }
}

template<class Container, int numPhases, int numComponents>
void testMSContainer()
{
    TestFluidSystem<numPhases, numComponents> fs;
    Container diffCoeff;

    std::cout << "Testing " << Dune::className(diffCoeff) << std::endl;

    // fill the container with coefficients
    diffCoeff.update([fs](int phaseIdx, int compIIdx, int compJIdx) {
        return fs.binaryDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
    });

    // we test that the internal mapping is unique and that the interface is symmetric
    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
    {
        for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx)
        {
            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            {
                if (compIIdx != compJIdx)
                {
                    const auto refD = fs.binaryDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
                    const auto D = diffCoeff(phaseIdx, compIIdx, compJIdx);
                    if (D != refD)
                        DUNE_THROW(Dune::Exception,
                            "Retrieved incorrect diffusion coefficient for (p,c0,c1) = ("
                            << phaseIdx << "," << compIIdx << "," << compJIdx << "). "
                            << "Result: " << D << ", expected result: " << refD);
                }
            }
        }
    }
}

template<int numComponents, int numPhases>
using FickDC = Dumux::FickianDiffusionCoefficients<int, numComponents, numPhases>;

template<int numComponents, int numPhases>
using MSDC = Dumux::MaxwellStefanDiffusionCoefficients<int, numComponents, numPhases>;

int main(int argc, char* argv[])
{
    testFickianContainer<FickDC<1, 5>, 1, 5>();
    testFickianContainer<FickDC<2, 2>, 2, 2>();
    testFickianContainer<FickDC<2, 3>, 2, 3>();
    testFickianContainer<FickDC<3, 3>, 3, 3>();
    testFickianContainer<FickDC<3, 8>, 3, 8>();
    testFickianContainer<FickDC<3, 2>, 3, 2>();

    testMSContainer<MSDC<1, 1>, 1, 1>();
    testMSContainer<MSDC<1, 2>, 1, 2>();
    testMSContainer<MSDC<1, 3>, 1, 3>();
    testMSContainer<MSDC<2, 2>, 2, 2>();
    testMSContainer<MSDC<3, 8>, 3, 8>();
    testMSContainer<MSDC<3, 2>, 3, 2>();

    return 0;
}
