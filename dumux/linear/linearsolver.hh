/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Base class for all linear solvers
 */
#ifndef DUMUX_LINEAR_SOLVER_HH
#define DUMUX_LINEAR_SOLVER_HH

namespace Dumux {

/*!
 * \brief This is the base class for all linear solvers available in
 *        Dumux
 */
template <class Matrix, class Vector>
class LinearSolver
{
public:
    LinearSolver()
    {
    };

    /*!
     * \brief Prepare to solve a linear system of equations.
     * 
     * This method allocates space an does the necessarry
     * communication before actually calling the solve() method.  As
     * long as the structure of the linear system does not change, the
     * solve method can be called arbitrarily often.
     */
    void prepare(const Matrix &M, Vector &x, const Vector &b)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "LinearSolver::prepare() not overloaded by the derived class");
    };

    /*!
     * \brief Actually solve the linear system of equations. 
     *
     * \return true if the residual reduction could be achieved, else false.
     */
    bool solve(const Matrix &M, Vector &x, const Vector &b, double residReduction)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "LinearSolver::solve() not overloaded by the derived class");
    };

    /*!
     * \brief Clean up after the last call of the solve() method
     */
    void cleanup()
    {
        DUNE_THROW(Dune::NotImplemented,
                   "LinearSolver::prepare() not overloaded by the derived class");
    };
};

} // namespace Dumux

#endif
