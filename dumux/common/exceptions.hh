// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Some exceptions thrown in DuMu<sup>x</sup>
 */
#ifndef DUMUX_EXCEPTIONS_HH
#define DUMUX_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

#include <string>

namespace Dumux {
/*!
 * \ingroup Core
 * \brief Exception thrown if a fixable numerical problem occurs.
 *
 * (e.g. time step too big, etc.)
 */
class NumericalProblem : public Dune::Exception
{
public:
    // copy constructor
    NumericalProblem(const NumericalProblem &v)
        : Dune::Exception(v)
    {}

    // default constructor
    NumericalProblem()
    {}

    // constructor with error message
    explicit NumericalProblem(const std::string &s)
    { this->message(s); }
};

/*!
 * \ingroup Core
 * \brief Exception thrown if a run-time parameter is not specified correctly.
 */
class ParameterException : public Dune::Exception
{
public:
    // copy constructor
    ParameterException(const ParameterException &v)
        : Dune::Exception(v)
    {}

    // default constructor
    ParameterException()
    {}

    // constructor with error message
    explicit ParameterException(const std::string &s)
    { this->message(s); }
};

} // namespace Dumux

#endif
