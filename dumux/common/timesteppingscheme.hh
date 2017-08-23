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
 * \brief Manages the handling of time dependent problems
 */
#ifndef DUMUX_TIME_STEPPING_SCHEME_HH
#define DUMUX_TIME_STEPPING_SCHEME_HH

#include "parameters.hh"

namespace Dumux
{

/*!
 * \brief The abstract time stepping parameter interface
 */
template <class Scalar>
class TimeSteppingParams
{

public:
    //! Pure abstract base classes have virtual destructor
    virtual ~TimeSteppingScheme () {}

    //! Returns if the time stepping scheme is implicit
    virtual constexpr bool implicit() const = 0;

    //! The number of stages in the time stepping scheme
    virtual std::size_t numStages() const = 0;

    //! The a parameters of the time stepping scheme
    virtual Scalar a(int stage, int i) const = 0;

    //! The b parameters of the time stepping scheme
    virtual Scalar b(int stage, int i) const = 0;

    //! The d parameters of the time stepping scheme
    virtual Scalar d(int i) const = 0;

    //! The name of the time stepping scheme
    virtual std::string name() const = 0;
};

/*!
 * \brief The time stepping scheme class
 */
template <class Scalar>
class TimeSteppingScheme
{

public:
    TimeSteppingScheme(std::shared_ptr<TimeSteppingParams> method)
    : method_(method)
    {}

private:
    std::shared_ptr<TimeSteppingParams> method_;

};

} // end namespace Dumux

#endif
