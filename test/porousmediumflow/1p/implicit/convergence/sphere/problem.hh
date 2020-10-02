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
 * \ingroup OnePTests
 * \brief The problem setup for the convergence test with analytic solution
 */
#ifndef DUMUX_CONVERGENCE_SPHERE_TEST_ONEP_PROBLEM_HH
#define DUMUX_CONVERGENCE_SPHERE_TEST_ONEP_PROBLEM_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The problem setup for the convergence test with analytic solution
 */
template <class TypeTag>
class ConvergenceSphereProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    ConvergenceSphereProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos)[0]; }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    Dune::FieldVector<Scalar, 1> analyticalSolution(const GlobalPosition& globalPos) const
    {
        // fundamental solution of laplace equation in 3d
        const auto r = globalPos.two_norm();
        return { 1.0/(4.0*M_PI*r) };
    }

    Dune::FieldVector<Scalar, 3> analyticalGradient(const GlobalPosition& globalPos) const
    {
        // fundamental solution of laplace equation in 3d
        const auto r2 = globalPos.two_norm2();
        const auto d = -1.0/(4.0*M_PI*r2);
        return { d, d, d };
    }

    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
};

} // end namespace Dumux

#endif
