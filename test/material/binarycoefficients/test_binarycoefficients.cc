// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief This test makes sure that the programming interface is
 *        observed by all binary coefficients.
 */

#include <config.h>

#include <dumux/material/binarycoefficients/air_mesitylene.hh>
#include <dumux/material/binarycoefficients/air_xylene.hh>
#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/binarycoefficients/h2o_ch4.hh>
#include <dumux/material/binarycoefficients/h2o_constant.hh>
#include <dumux/material/binarycoefficients/h2o_heavyoil.hh>
#include <dumux/material/binarycoefficients/h2o_mesitylene.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>
#include <dumux/material/binarycoefficients/h2o_co2.hh>
#include <dumux/material/binarycoefficients/h2o_ar.hh>
#include <dumux/material/binarycoefficients/h2o_xylene.hh>
#include <dumux/material/binarycoefficients/n2_o2.hh>

#include <dumux/material/components/co2.hh>
#include <dumux/material/components/defaultco2table.hh>

template<class Scalar, class BinaryCoefficients>
int checkBinaryCoefficients()
{
    // Define a reference temperature and pressure for evaluation:
    const Scalar t = 293.15;
    const Scalar p = 1.013e5;

    const std::string className = Dune::className<BinaryCoefficients>();
    std::string collectedErrors;

    // Test the function gasDiffCoeff:
    try
    {
        [[maybe_unused]] Scalar val = BinaryCoefficients::gasDiffCoeff(t, p);
    } catch (Dune::NotImplemented&)
    {
        std::cout << "warning: " << className << "::gasDiffCoeff() is not implemented." << std::endl;
    } catch (const std::exception& e)
    {
        collectedErrors += "error: " + className + "::gasDiffCoeff() throws exception: "
                           + std::string(e.what()) + "\n";
    }

    // Test the function liquidDiffCoeff:
    try
    {
        [[maybe_unused]] Scalar val = BinaryCoefficients::liquidDiffCoeff(t, p);
    } catch (Dune::NotImplemented&)
    {
        std::cout << "warning: " << className << "::liquidDiffCoeff() is not implemented." << std::endl;
    } catch (const std::exception& e)
    {
        collectedErrors += "error: " + className + "::liquidDiffCoeff() throws exception: "
                           + std::string(e.what()) + "\n";
    }

    // Test the function henry:
    try
    {
        [[maybe_unused]] Scalar val = BinaryCoefficients::henry(t);
        std::cout << className << ": t=" << t << ", p=" << p << ", henry: " << val << std::endl;
    } catch (Dune::NotImplemented&)
    {
        std::cout << "warning: " << className << "::henry() is not implemented." << std::endl;
    } catch (const std::exception& e)
    {
        collectedErrors += "error: " + className + "::henry() throws exception: "
                           + std::string(e.what()) + "\n";
    }

    //Test the henryMixture of air
    try
    {
        [[maybe_unused]] Scalar val = BinaryCoefficients::henryMixture(t);
        std::cout << className << ": t=" << t << ", p=" << p << ", henryMixture: " << val << std::endl;
    } catch (Dune::NotImplemented&)
    {
        std::cout << "warning: " << className << "::henryMixture() is not implemented." << std::endl;
    } catch (const std::exception& e)
    {
        collectedErrors += "error: " + className + "::henryMixture() throws exception: "
                           + std::string(e.what()) + "\n";
    }



    std::cout << collectedErrors;
    if (collectedErrors.empty()) // success
        return 0;
    else
        return 1;
}

int main()
{
    using namespace Dumux;
    using namespace Dumux::BinaryCoeff;
    using Scalar = double;

    int success = 0;

//    success += checkBinaryCoefficients<Scalar, Air_Mesitylene>();
//    success += checkBinaryCoefficients<Scalar, Air_Xylene>();
//    success += checkBinaryCoefficients<Scalar, Brine_CO2<Scalar, Components::CO2<Scalar, GeneratedCO2Tables::CO2Tables>>>();
   success += checkBinaryCoefficients<Scalar, H2O_Air>();
//    success += checkBinaryCoefficients<Scalar, H2O_CH4>();
//    success += checkBinaryCoefficients<Scalar, H2O_Component<Scalar, Components::Constant<0, Scalar>>>();
//    success += checkBinaryCoefficients<Scalar, H2O_HeavyOil>();
//    success += checkBinaryCoefficients<Scalar, H2O_Mesitylene>();
   success += checkBinaryCoefficients<Scalar, H2O_N2>();
   success += checkBinaryCoefficients<Scalar, H2O_O2>();
   success += checkBinaryCoefficients<Scalar, H2O_CO2>();
   success += checkBinaryCoefficients<Scalar, H2O_AR>();
//    success += checkBinaryCoefficients<Scalar, H2O_Xylene>();
//    success += checkBinaryCoefficients<Scalar, N2_O2>();

    return success;
}
