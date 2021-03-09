#include <config.h>

#include <cmath>
#include <iterator>
#include <string>
#include <string_view>

#include <dune/common/exceptions.hh>

#include <dumux/io/format.hh>

void testString(const std::string_view result, const std::string_view expected)
{
    if (result != expected)
        DUNE_THROW(Dune::Exception, "Unexpected result: " << result << ", expected " << expected);
}

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // format
    // format string is like in Python and documented here:
    // https://en.cppreference.com/w/cpp/utility/format/formatter#Standard_format_specification
    {
        testString(Fmt::format("Number {}", 42), "Number 42");
        testString(Fmt::format("Numbers {0}, {0}, {1}", 42, 45), "Numbers 42, 42, 45");
        testString(Fmt::format("Numbers {0}, {0}, {{{1}}}", 42, 45), "Numbers 42, 42, {45}");
        testString(Fmt::format("{:.3f}", M_PI), "3.142");
        testString(Fmt::format("{:.2e}", M_PI), "3.14e+00");
        testString(Fmt::format("{:6}", "x"), "x     ");
    }

    // format_to
    {
        std::string buffer;
        Fmt::format_to(std::back_inserter(buffer), "{:.2f}!", M_PI);
        testString(buffer, "3.14!");
    }

    // format_to_n (from https://en.cppreference.com/w/cpp/utility/format/format_to_n)
    {
        char buffer[64];
        const auto result = Fmt::format_to_n(
            buffer, std::size(buffer),
            "Hubble's H{0} {1} {2} miles/sec/mpc.",
            "\u2080", "\u2245", 42
        );

        testString(std::string_view{buffer, result.size}, "Hubble's H₀ ≅ 42 miles/sec/mpc.");
    }

    // formatted_size
    {
        std::string_view fmt{"{}"};
        int num = 42;
        const auto bufferSize = Fmt::formatted_size(fmt, num);
        if (bufferSize != 2)
            DUNE_THROW(Dune::Exception, "Unexpected computed buffer size: " << bufferSize << ", expected 2 for storing '42'");
    }


    return 0;
}
