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
 * \ingroup Common
 * \brief A high-level interface for solvers of PDEs.
 */
#ifndef DUMUX_COMMON_PDE_SOLVER_HH
#define DUMUX_COMMON_PDE_SOLVER_HH

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A high-level interface for solvers of PDEs.
 * \tparam V The variables of the system of equations
 */
template<typename V>
class PDESolver
{
public:
    using Variables = V;

    virtual ~PDESolver() = default;

    //! Solve the system and store the result in the given variables object.
    bool solve(Variables& vars)
    { return solve_(vars); }

private:
    virtual bool solve_(Variables&) = 0;
};

} // namespace Dumux

#endif
