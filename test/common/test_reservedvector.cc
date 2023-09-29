//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <dune/common/exceptions.hh>
#include <dumux/common/reservedvector.hh>

static constexpr int numValues = 1000;

template<std::size_t N>
void testValues(const Dumux::ReservedVector<int, N>& v)
{
    for (int i = 0; i < numValues; i++)
        if (v[i] != i)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
}

template<std::size_t N>
void test(Dumux::ReservedVector<int, N>&& v)
{
    v.clear();
    for (int i = 0; i < numValues; i++)
        v.push_back(i);
    testValues(v);
    auto cpy = v; testValues(cpy);
    auto moved = std::move(v); testValues(cpy);
    Dumux::ReservedVector moved_cted{std::move(moved)}; testValues(moved_cted);
}

int main(int argc, char* argv[])
{
    test(Dumux::ReservedVector<int, 2000>{});
    test(Dumux::ReservedVector<int, 1000>{});
    test(Dumux::ReservedVector<int, 10>{});

    return 0;
}
