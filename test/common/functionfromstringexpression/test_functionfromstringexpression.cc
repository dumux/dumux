//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/functionfromstringexpression.hh>

int main(int argc, char* argv[])
{
    using namespace Dumux;

    constexpr std::size_t numSamples = 100;
    const auto values0 = linspace(-M_PI, M_PI, numSamples);
    const auto values1 = linspace(-1.0, 1.0, numSamples);

    const auto test1Arg = [&](auto&& f1, auto&& f2)
    {
        for (int i = 0; i < numSamples; ++i)
            if (Dune::FloatCmp::ne<double>(f1(values0[i]), f2(values0[i])))
                DUNE_THROW(Dune::Exception, "Results do not match: "
                    << f1(values0[i]) << " " << f2(values0[i]));

    };

    const auto test2Args = [&](auto&& f1, auto&& f2)
    {
        for (int i = 0; i < numSamples; ++i)
            if (Dune::FloatCmp::ne<double>(f1(values0[i], values1[i]), f2(values0[i], values1[i])))
                DUNE_THROW(Dune::Exception, "Results do not match: "
                    << f1(values0[i], values1[i]) << " " << f2(values0[i], values1[i]));
    };

    // function with one argument
    {
        const std::string funcStr = "sin(t)";
        FunctionFromStringExpression<1> func(funcStr, "t");
        test1Arg(func, [](auto t) { return std::sin(t); });
    }

    // function with two arguments
    {
        const std::string funcStr = "x*sin(t)";
        FunctionFromStringExpression<2> func(funcStr, "xt");
        test2Args(func, [](auto x, auto t) { return x*std::sin(t); });
    }

    // function with longer variable names
    {
        const std::string funcStr = "foo*sin(bar)";
        FunctionFromStringExpression<2> func(funcStr, std::array<std::string, 2>{{"foo", "bar"}});
        test2Args(func, [](auto foo, auto bar) { return foo*std::sin(bar); });
    }

    // function with expression-local variables
    {
        const std::string funcStr = "var b:=2; x*sin(t) + b + 24;";
        FunctionFromStringExpression<2> func(funcStr, "xt");
        test2Args(func, [](auto x, auto t) { double b = 2; return x*std::sin(t) + b + 24; });
    }

    // faulty expression that has to throw an exception
    {
        const std::string funcStr = "foo*sin(bar)*bla";
        bool error = false;
        try { FunctionFromStringExpression<2> func(funcStr, std::array<std::string, 2>{{"foo" "bar"}}); }
        catch (Dune::IOError& e) { error = true; }
        if (!error)
            DUNE_THROW(Dune::Exception, "Faulty expression did not throw: " << funcStr);
    }

    return 0;
}
