// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Writes a container to file
 * \param v The container, requires begin() and end() method
 * \param filename The filename to write to
 * \param floatPrecision The total number of digits stored, including decimal
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
 * \brief Read an input stream into a container
 * \param stream A standard input stream
 * \tparam Container The container type requires begin(), end(), push_back() functions
 *                   and Container::value_type requires operator>>.
 */
template<typename Container>
Container readStreamToContainer(std::istream& stream)
{
    Container v;
    std::istream_iterator<typename Container::value_type> it(stream);
    std::copy(it, std::istream_iterator<typename Container::value_type>(), std::back_inserter(v));
    return v;
}

/*!
 * \brief Read a simple text file into a container
 * \param filename The filename to write to
 * \tparam Container The container type requires begin(), end(), push_back() functions
 *                   and Container::value_type requires operator>>.
 *
 * usage: auto v = readFileToContainer<std::vector<double>>("myvector.txt");
 */
template<typename Container>
Container readFileToContainer(const std::string& filename)
{
    std::ifstream infile(filename, std::ios::in);
    if (!infile)
        DUNE_THROW(Dune::IOError, "Could not open file: " << filename);
    return readStreamToContainer<Container>(infile);
}

} // end namespace Dumux

#endif
