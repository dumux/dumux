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
 * \brief Test writing and reading sequence container to and from file
 */
#include <config.h>

#include <algorithm>
#include <vector>
#include <list>
#include <deque>
#include <array>
#include <initializer_list>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dumux/io/container.hh>

namespace Dumux {

template<typename T>
bool testContainerIO(const T& c0, int floatPrecision = 6)
{
    writeContainerToFile(c0, "container.txt", floatPrecision);
    auto c1 = readFileToContainer<T>("container.txt");
    return std::equal(c0.begin(), c0.end(), c1.begin());
}

template<typename T>
bool testContainerIO(const std::initializer_list<typename T::value_type>& init)
{
    T c0(init);
    writeContainerToFile(c0, "container.txt");
    auto c1 = readFileToContainer<T>("container.txt");
    return std::equal(c0.begin(), c0.end(), c1.begin());
}

template<typename T>
bool testContainerWriter(const T& c0)
{
    writeContainerToFile(c0, "container.txt");
    auto c1 = readFileToContainer<std::vector<std::decay_t<decltype(*c0.begin())>>>("container.txt");
    return std::equal(c0.begin(), c0.end(), c1.begin());
}

} // end namespace Dumux

////////////////////////
// the main function
////////////////////////
int main()
{
    bool passed = true;

    // we can read and write into
    // std::vector, std::list, std::deque
    // we can even read and write containers of FieldVectors (very convenient!)
    auto doublei = {5.3, 6.1, 7.2, 5.66, 2.89, 8.123};
    passed = passed && Dumux::testContainerIO<std::vector<double>>(doublei);
    passed = passed && Dumux::testContainerIO<std::list<double>>(doublei);
    passed = passed && Dumux::testContainerIO<std::deque<double>>(doublei);
    if (!passed) return 1;

    auto doublepreci = {1.23456789123456, 1.23456789123456, 9.87654321987654};
    passed = passed && Dumux::testContainerIO<std::vector<double>>(doublepreci, 15);
    passed = passed && !Dumux::testContainerIO<std::vector<double>>(doublepreci, 7);
    if (!passed) return 1;

    auto inti = {5, 6, 7, 5, 2, 8};
    passed = passed && Dumux::testContainerIO<std::vector<int>>(inti);
    passed = passed && Dumux::testContainerIO<std::list<int>>(inti);
    passed = passed && Dumux::testContainerIO<std::deque<int>>(inti);
    if (!passed) return 1;

    std::initializer_list<std::string> stringi = {"5", "6", "7", "5", "2", "8"};
    passed = passed && Dumux::testContainerIO<std::vector<std::string>>(stringi);
    passed = passed && Dumux::testContainerIO<std::list<std::string>>(stringi);
    passed = passed && Dumux::testContainerIO<std::deque<std::string>>(stringi);
    if (!passed) return 1;

    auto fvector2i = {Dune::FieldVector<double, 2>(0.0), Dune::FieldVector<double, 2>(1.0), Dune::FieldVector<double, 2>(2.0)};
    passed = passed && Dumux::testContainerIO<std::vector<Dune::FieldVector<double, 2>>>(fvector2i);
    passed = passed && Dumux::testContainerIO<std::list<Dune::FieldVector<double, 2>>>(fvector2i);
    passed = passed && Dumux::testContainerIO<std::deque<Dune::FieldVector<double, 2>>>(fvector2i);
    if (!passed) return 1;

    auto fvector3i = {Dune::FieldVector<double, 3>(0.0), Dune::FieldVector<double, 3>(1.0), Dune::FieldVector<double, 3>(2.0)};
    passed = passed && Dumux::testContainerIO<std::vector<Dune::FieldVector<double, 3>>>(fvector3i);
    passed = passed && Dumux::testContainerIO<std::list<Dune::FieldVector<double, 3>>>(fvector3i);
    passed = passed && Dumux::testContainerIO<std::deque<Dune::FieldVector<double, 3>>>(fvector3i);
    if (!passed) return 1;

    // we can write also std::arrays (all container providing begin and end)
    passed = passed && Dumux::testContainerWriter<std::array<double, 2>>({{1.0, 2.0}});
    passed = passed && Dumux::testContainerWriter<std::array<int, 2>>({{1, 2}});
    passed = passed && Dumux::testContainerWriter<std::array<std::string, 2>>({{"1.0", "2.0"}});
    passed = passed && Dumux::testContainerWriter<std::array<Dune::FieldVector<double, 3>, 2>>(
                        std::array<Dune::FieldVector<double, 3>, 2>{{Dune::FieldVector<double, 3>(0.0), Dune::FieldVector<double, 3>(1.0)}});
    if (!passed) return 1;

    return 0;
}
