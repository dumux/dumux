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
int checkDiffusionCoefficients()
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

    std::cout << collectedErrors;
    if (collectedErrors.empty()) // success
        return 0;
    else
        return 1;
}

struct VoidTag{};

template<class Scalar, class BinaryCoefficients, class Implementation=VoidTag>
int checkHenryCoefficients(Implementation impl = Implementation{})
{
    // Define a reference temperature and pressure for evaluation:
    const Scalar t = 293.15;
    const Scalar p = 1.013e5;

    const std::string className = Dune::className<BinaryCoefficients>();
    std::string collectedErrors;

    try
    {
        if constexpr (std::is_same_v<Implementation, VoidTag>) {
            [[maybe_unused]] Scalar val = BinaryCoefficients::henry(t);
            std::cout << className << ": t=" << t << ", p=" << p << ", henry: " << val << std::endl;
        } else {
            [[maybe_unused]] Scalar val = BinaryCoefficients::henry(t, impl);
            std::cout << className << ": t=" << t << ", p=" << p << ", henry with specific implementation: " << val << std::endl;
        }
    } catch (Dune::NotImplemented&)
    {
        std::cout << "warning: " << className << "::henry() is not implemented." << std::endl;
    } catch (const std::exception& e)
    {
        collectedErrors += "error: " + className + "::henry() throws exception: "
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

    success += checkDiffusionCoefficients<Scalar, Air_Mesitylene>();
    success += checkDiffusionCoefficients<Scalar, Air_Xylene>();
    success += checkDiffusionCoefficients<Scalar, Brine_CO2<Scalar, Components::CO2<Scalar, GeneratedCO2Tables::CO2Tables>>>();
    success += checkDiffusionCoefficients<Scalar, H2O_Air>();
    success += checkDiffusionCoefficients<Scalar, H2O_CH4>();
    success += checkDiffusionCoefficients<Scalar, H2O_Component<Scalar, Components::Constant<0, Scalar>>>();
    success += checkDiffusionCoefficients<Scalar, H2O_HeavyOil>();
    success += checkDiffusionCoefficients<Scalar, H2O_Mesitylene>();
    success += checkDiffusionCoefficients<Scalar, H2O_N2>();
    success += checkDiffusionCoefficients<Scalar, H2O_O2>();
    success += checkDiffusionCoefficients<Scalar, H2O_CO2>();
    success += checkDiffusionCoefficients<Scalar, H2O_AR>();
    success += checkDiffusionCoefficients<Scalar, H2O_Xylene>();
    success += checkDiffusionCoefficients<Scalar, N2_O2>();

    success += checkHenryCoefficients<Scalar, Air_Mesitylene>();
    success += checkHenryCoefficients<Scalar, Air_Xylene>();
    success += checkHenryCoefficients<Scalar, H2O_Air>();
    success += checkHenryCoefficients<Scalar, H2O_Air>(BinaryCoeff::Detail::H2O_Air::FinsterleImplementation{});
    success += checkHenryCoefficients<Scalar, H2O_Air>(BinaryCoeff::Detail::H2O_Air::IdealMixtureImplementation{});
    success += checkHenryCoefficients<Scalar, H2O_Air>(BinaryCoeff::Detail::H2O_Air::TchobanoglousCubicImplementation{});
    success += checkHenryCoefficients<Scalar, H2O_CH4>();
    success += checkHenryCoefficients<Scalar, H2O_Mesitylene>();
    success += checkHenryCoefficients<Scalar, H2O_N2>();
    success += checkHenryCoefficients<Scalar, H2O_O2>();
    success += checkHenryCoefficients<Scalar, H2O_CO2>();
    success += checkHenryCoefficients<Scalar, H2O_AR>();
    success += checkHenryCoefficients<Scalar, H2O_Xylene>();
    success += checkHenryCoefficients<Scalar, N2_O2>();

    return success;
}
