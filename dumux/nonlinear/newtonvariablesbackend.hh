// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup Nonlinear
 * \brief Backends for the Newton for different variable types
 */
#ifndef DUMUX_NEWTON_VARIABLES_BACKEND_HH
#define DUMUX_NEWTON_VARIABLES_BACKEND_HH

#include <type_traits>
#include <dune/common/typetraits.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup Nonlinear
 * \brief Class providing operations with variables needed in the Newton solver
 */
template<class Variables, class Enable = Variables>
class NewtonVariablesBackend;

/*!
 * \file
 * \ingroup Nonlinear
 * \brief Class providing Newton operations for scalar/number types
 */
template<class Scalar>
class NewtonVariablesBackend<Scalar, std::enable_if_t<Dune::IsNumber<Scalar>::value, Scalar>>
{
public:
    using Variables = Scalar; //!< the type of the variables object
    using DofVector = Scalar; //!< the type of the dofs parametrizing the variables object

    static std::size_t size(const DofVector& d)
    { return 1; }

    static DofVector makeDofVector(const DofVector& d)
    { return d; }

    static DofVector makeZeroDofVector(std::size_t size)
    { return 0.0; }

    static Scalar makeDofVectorForSolver(const DofVector& d)
    { return d; }

    static Scalar makeZeroDofVectorForSolver(std::size_t size)
    { return 0.0; }

    static Scalar reconstructDofVectorFromSolver(const Scalar& d)
    { return d; }

    static Scalar maxRelativeShift(const DofVector& previous, const DofVector& delta)
    {
        const auto current = previous - delta;
        using std::abs; using std::max;
        auto error = abs(previous - current);
        error /= max<Scalar>(1.0, abs(previous + current)/2.0);
        return error;
    }

    // //! operations on variables
    // static DofVector& getDofVector(Variables& v)
    // { return v; }
};

} // end namespace Dumux

#endif
