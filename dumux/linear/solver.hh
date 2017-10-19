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
 * \brief Dumux linear solver backend
 */
#ifndef DUMUX_SOLVER_BACKEND_HH
#define DUMUX_SOLVER_BACKEND_HH

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ImplicitIsBox);
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(ElementMapper);
}

/*!
 * \ingroup Linear
 * \brief Base class for linear solvers
 */
template<class TypeTag>
class LinearSolver
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using DofMapper = std::conditional_t<GET_PROP_VALUE(TypeTag, ImplicitIsBox),
                                         typename GET_PROP_TYPE(TypeTag, VertexMapper),
                                         typename GET_PROP_TYPE(TypeTag, ElementMapper)>;

public:
    LinearSolver (const GridView& gridView, const DofMapper& mapper)
    : gridView_(gridView)
    , mapper_(mapper)
    {}

    //! solve the linear system Ax = b
    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    { return false; };

    //! the name of the linear solver
    std::string name() const
    { return "undefined solver"; }

protected:
    //! the grid view for communication in parallel
    const GridView& gridView() const
    { return gridView_; }

    //! the dof mapper for communication in parallel
    const DofMapper& dofMapper() const
    { return mapper_; }

private:
    const GridView gridView_;
    const DofMapper mapper_;
};

} // end namespace Dumux

#endif
