#include <config.h>

#include <iostream>
#include <stdexcept>

#include <string>
#include <dumux/common/stringutilities.hh>

int main(int argc, char* argv[]) try
{
    const auto checkTokens = [](const std::string& name, const auto& tokens, const std::vector<std::string>& ref)
    {
        if (tokens.size() != ref.size())
            throw std::runtime_error(name + "Wrong number of tokens!");

        for (int i = 0; i < ref.size(); ++i)
            if (tokens[i] != ref[i])
                throw std::runtime_error(name + "Wrong token! Expected: '" + ref[i] + "' got '" + std::string(tokens[i]) + "'");

        std::cout << name << "passed" << std::endl;
    };

    std::string str;

    // tokenzie
    std::cout << "------------\nDumux::tokenize\n---------------\n";

    str = "bla&foo&bar";
    checkTokens("Test 1: ", Dumux::tokenize(str, "&"), {"bla", "foo", "bar"});

    str = "&1&";
    checkTokens("Test 2: ", Dumux::tokenize(str, "&"), {"1"});

    str = "&&&1&";
    checkTokens("Test 3: ", Dumux::tokenize(str, "&"), {"1"});

    str = "&&1&2&&3&";
    checkTokens("Test 4: ", Dumux::tokenize(str, "&"), {"1", "2", "3"});

    str = "    a      b     ";
    checkTokens("Test 5: ", Dumux::tokenize(str, " "), {"a", "b"});

    str = "    a      b     c";
    checkTokens("Test 6: ", Dumux::tokenize(str, " "), {"a", "b", "c"});

    str = "ababababab";
    checkTokens("Test 7: ", Dumux::tokenize(str, "ab"), {});

    str = "| hello | world |\n\n";
    checkTokens("Test 8: ", Dumux::tokenize(str, "| \n"), {"hello", "world"});


    // split
    std::cout << "------------\nDumux::split\n---------------\n";

    str = "bla&foo&bar";
    checkTokens("Test 1: ", Dumux::split(str, "&"), {"bla", "foo", "bar"});

    str = "fooborfoo";
    checkTokens("Test 2: ", Dumux::split(str, "bor"), {"foo", "foo"});

    str = "fooborfoo";
    checkTokens("Test 3a: ", Dumux::split(str, "foo"), {"", "bor", ""});
    checkTokens("Test 3b: ", Dumux::split(str, "foo", true), {"bor"});

    return 0;
}
catch (const std::exception& e)
{
    std::cout << e.what();
    return 1;
}
