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
 * \ingroup InputOutput
 * \brief Free functions to write and read a sequence container to and from a file
 * \note Reading should work for all sequence containers providing begin, end, and push_back
 *       (e.g. std::vector, std::deque, std::list), so not for e.g. std::array.
 *       Writing only needs begin and end member functions returning iterators.
 */
#ifndef DUMUX_IO_CONTAINER_HH
#define DUMUX_IO_CONTAINER_HH

#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream>
#include <iterator>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Writes a container to file
 * \param v The container, requires begin() and end() method
 * \param filename The filename to write to
 *
 * usage: std::vector<double> v(5, 0.0); writeContainerToFile(v, "myvector.txt");
 */
template<typename Container>
void writeContainerToFile(const Container& v,
                          const std::string& filename,
                          int floatPrecision = 6)
{
    std::ofstream outfile(filename, std::ios::out);
    outfile << std::scientific << std::setprecision(floatPrecision);
    std::ostream_iterator<typename Container::value_type> it(outfile, "\n");
    std::copy(v.begin(),v.end(), it);
}

/*!
 * \brief Read a simple text file into a container
 * \param filename The filename to write to
 * \tparam Container  The container type, requires begin(), end(), push_back() method
 *
 * usage: auto v = readFileToContainer<std::vector<double>>("myvector.txt");
 */
template<typename Container>
Container readFileToContainer(const std::string& filename)
{
    Container v;
    std::ifstream infile(filename, std::ios::in);
    std::istream_iterator<typename Container::value_type> it(infile);
    std::copy(it, std::istream_iterator<typename Container::value_type>(), std::back_inserter(v));
    return v;
}

} // end namespace Dumux

#endif
