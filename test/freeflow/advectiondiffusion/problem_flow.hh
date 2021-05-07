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
 * \ingroup AdvectionDiffusionTests
 * \brief Definition of a problem for the advection diffusion transport problem:
 */
#ifndef DUMUX_ADVECTION_DIFFUSION_FLOW_TEST_PROBLEM_HH
#define DUMUX_ADVECTION_DIFFUSION_FLOW_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * Flow from left to right in a two-dimensional channel is considered. At the inlet (left),
 * fixed values for velocity are set, while at the outlet (right), a fixed pressure
 * boundary condition is used. The channel is confined by solid walls at the top and bottom
 * of the domain which corresponds to no-slip/no-flow conditions.
 * For the non-isothermal test, water of increased temperature is injected at the inlet
 * while the walls are fully isolating.
 */
template <class TypeTag>
class FlowProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr auto dimWorld = GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

public:
    FlowProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { inletVelocity_ = getParam<Scalar>("Problem.InletVelocity"); }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if (isInlet_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else if (isOutlet_(globalPos))
            values.setDirichlet(Indices::pressureIdx);
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        return values;
    }

    /*!
      * \brief Evaluates the boundary conditions for a Dirichlet control volume.
      *
      * \param globalPos The center of the finite volume which ought to be set.
      */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        if (isInlet_(globalPos))
            values[Indices::velocityXIdx] = parabolicProfile(globalPos[1], inletVelocity_);
        return values;
    }

    /*!
     * \brief A parabolic velocity profile.
     *
     * \param y The position where the velocity is evaluated.
     * \param vMax The profile's maxmium velocity.
     */
    Scalar parabolicProfile(const Scalar y, const Scalar vMax) const
    {
        const Scalar yMin = this->gridGeometry().bBoxMin()[1];
        const Scalar yMax = this->gridGeometry().bBoxMax()[1];
        return  vMax * (y - yMin)*(yMax - y) / (0.25*(yMax - yMin)*(yMax - yMin));
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        values[Indices::pressureIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;
        values[Indices::velocityXIdx] = 0.0;

        return values;
    }

private:

    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    static constexpr Scalar eps_=1e-6;
    Scalar inletVelocity_;

};
} // end namespace Dumux

#endif
