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
/**
 * \file
 * \ingroup OnePTests
 * \brief Test for the OnePModel in combination with the NI model for a simple problem.
 *
 * The simulation domain is the unit square. While a pressure gradient from left
 * to right is imposed by corresponding Dirichlet conditions, the prescribed
 * temperature is the same on both sides. Therefore, a uniform constant temperature
 * distribution should be expected.
 */

#ifndef DUMUX_TEST_1PNI_SIMPLE_PROBLEM_HH
#define DUMUX_TEST_1PNI_SIMPLE_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Test for the OnePModel in combination with the NI model
 */
template <class TypeTag>
class OnePNISimpleProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

public:
    OnePNISimpleProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        rightPressure_ = getParam<Scalar>("Problem.RightPressure", 1e5);
        injectionRate_ = getParam<Scalar>("Problem.InjectionRate", 0.0);
        multiplier_ = getParam<Scalar>("Problem.Multiplier", 1.0);
    }

    const std::string& name() const
    { return name_; }


    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            return {1e5, 300.0};
        else
            return {rightPressure_, 300.0};
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        if (globalPos[0] < 0.6 && globalPos[0] > 0.4 && globalPos[1] < 0.6 && globalPos[1] > 0.4)
        {
            values[Indices::conti0EqIdx] = injectionRate_;
            values[Indices::conti0EqIdx] = injectionRate_*multiplier_;
        }

        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return dirichletAtPos(globalPos); }

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    Scalar rightPressure_;
    Scalar injectionRate_;
    Scalar multiplier_;
};

} // end namespace Dumux

#endif
