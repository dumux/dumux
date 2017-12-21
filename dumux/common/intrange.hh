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
 * \ingroup Common
 * \brief A class to create an iterable integer range
 */
#ifndef DUMUX_INTEGER_RANGE_HH
#define DUMUX_INTEGER_RANGE_HH

#include <cassert>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief This class generates a IntRange [a,b) which can be used in a for loop, e.g.
 *        for(auto i : IntRange(3) { ... i = 0, 1, 2 }   or
 *        for(auto i : IntRange(5, 8) { ... i = 5, 6, 7 }
 *        see: https://en.wikipedia.org/wiki/Generator_(computer_programming)
 */
class IntRange
{
public:
    // constructors
    IntRange(int end): last_(end), iter_(0) { assert(end > 0); }
    IntRange(int begin, int end): last_(end), iter_(begin) { assert(end > begin); }
    IntRange() = delete;

    // Iterable functions
     const IntRange& begin() const { return *this; }
     const IntRange& end() const { return *this; }

    // Iterator functions
    bool operator!=(const IntRange&) const { return iter_ < last_; }
    void operator++() { ++iter_; }
    int operator*() const { return iter_; }
private:
    int last_;
    int iter_;
};

} // namespace Dumux

#endif
