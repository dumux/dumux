// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief This is a program to test the flash calculation which uses
 *        non-linear complementarity problems (NCP)
 *
 * A flash calculation determines the pressures, saturations and
 * composition of all phases given the total mass (or, as in this case
 * the total number of moles) in a given amount of pore space.
 */
#include <config.h>
#include <functional>


#include <dune/common/timer.hh>
#include <dune/common/float_cmp.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>
#include <dumux/material/binarycoefficients/n2_o2.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/h2o.hh>


template <class Scalar, class GenericFunction, class OptimizedFunction>
static bool test(const std::string& name,
                 GenericFunction genericFunction,
                 OptimizedFunction optimizedFunction)
{
    const std::size_t size = 1e5;

    std::vector<Scalar> resultGeneric(size);
    std::vector<Scalar> resultOptimized(size);

    Dune::Timer timerGeneric;

    for(auto& x : resultGeneric)
        x = genericFunction();

    timerGeneric.stop();

    Dune::Timer timerOptimized;

    for(auto& x : resultOptimized)
        x = optimizedFunction();

    timerOptimized.stop();

    std::cout << name << " generic:   result " << resultGeneric[0] << ", time " << timerGeneric.elapsed() << " s" << std::endl;
    std::cout << name << " optimized: result " << resultOptimized[0] << ", time " << timerOptimized.elapsed() << " s" << std::endl;
    std::cout << std::endl;
    return Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(resultGeneric[0], resultOptimized[0], 1.0e-30);

}

int main()
{
    using namespace Dumux;
    using Scalar = double;
    const Scalar temperature = 300.0; // [K]
    const Scalar pressure = 1e5; // [Pa]

    std::vector<bool> testPassed;

    // lambdas for N2_O2
    auto genericGasDiffCoeff_N2_O2 = [temperature, pressure] ()
    {
        // molar masses [g/mol]
        const std::array<Scalar, 2> molarMass = {N2<Scalar>::molarMass()*1e3, O2<Scalar>::molarMass()*1e3};

        // atomic diffusion volumes
        const std::array<Scalar, 2> sigmaNu = {18.5 /* N2 */,  16.3 /* O2 */};

        return BinaryCoeff::fullerMethod(molarMass, sigmaNu, temperature, pressure);
    };

    auto optimizedGasDiffCoeff_N2_O2 = [temperature, pressure] ()
    { return BinaryCoeff::N2_O2::gasDiffCoeff(temperature, pressure); };



    // lambdas for H2O_O2
    auto genericGasDiffCoeff_H2O_O2 = [temperature, pressure] ()
    {
        // molar masses [g/mol]
        const std::array<Scalar, 2> molarMass = {H2O<Scalar>::molarMass()*1e3, O2<Scalar>::molarMass()*1e3};

        // atomic diffusion volumes
        const std::array<Scalar, 2> sigmaNu = {13.1 /* H2O */,  16.3 /* O2 */};

        return BinaryCoeff::fullerMethod(molarMass, sigmaNu, temperature, pressure);
    };

    auto optimizedGasDiffCoeff_H2O_O2  = [temperature, pressure] ()
    { return BinaryCoeff::H2O_O2::gasDiffCoeff(temperature, pressure); };

    // lambdas for H2O_N2
    auto genericGasDiffCoeff_H2O_N2 = [temperature, pressure] ()
    {
        // molar masses [g/mol]
        const std::array<Scalar, 2> molarMass = {H2O<Scalar>::molarMass()*1e3, N2<Scalar>::molarMass()*1e3};

        // atomic diffusion volumes
        const std::array<Scalar, 2> sigmaNu = {13.1 /* H2O */,  18.5 /* N2 */};

        return BinaryCoeff::fullerMethod(molarMass, sigmaNu, temperature, pressure);
    };

    auto optimizedGasDiffCoeff_H2O_N2  = [temperature, pressure] ()
    { return BinaryCoeff::H2O_N2::gasDiffCoeff(temperature, pressure); };

    testPassed.push_back(test<Scalar>("N2_O2", genericGasDiffCoeff_N2_O2, optimizedGasDiffCoeff_N2_O2));
    testPassed.push_back(test<Scalar>("H2O_O2", genericGasDiffCoeff_H2O_O2, optimizedGasDiffCoeff_H2O_O2));
    testPassed.push_back(test<Scalar>("H2O_N2", genericGasDiffCoeff_H2O_N2, optimizedGasDiffCoeff_H2O_N2));

    const int exitCode = std::all_of(testPassed.begin(), testPassed.end(), [](auto i) { return i; }) ? 0 : 1;


    return exitCode;
}
