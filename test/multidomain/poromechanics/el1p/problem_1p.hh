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
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the single-phase flow
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
#ifndef DUMUX_1P_SUB_PROBLEM_HH
#define DUMUX_1P_SUB_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The single-phase sub problem in the el1p coupled problem.
 */
template <class TypeTag>
class OnePSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // copy pressure index for convenience
    enum { pressureIdx = GetPropType<TypeTag, Properties::ModelTraits>::Indices::pressureIdx };

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    OnePSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                   const std::string& paramGroup = "OneP")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Returns the temperature within the domain in [K].
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    //! Evaluates the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0e5); }

    //! Evaluates source terms.
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar source = getParam<Scalar>("Problem.InjectionRate");
        if (globalPos[0] > 0.4 && globalPos[0] < 0.6 && globalPos[1] < 0.6 && globalPos[1] > 0.4)
            return NumEqVector(source);
        return NumEqVector(0.0);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

private:
    static constexpr Scalar eps_ = 1.0e-6;
    std::string problemName_;
};

} // end namespace Dumux

#endif
