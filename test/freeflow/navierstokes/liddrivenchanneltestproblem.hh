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
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model
 */
#ifndef DUMUX_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

namespace Dumux
{
template <class TypeTag>
class ChannelTestProblem;

namespace Properties
{
#if !NONISOTHERMAL
NEW_TYPE_TAG(ChannelTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokes));
#else
NEW_TYPE_TAG(ChannelTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNI));
#endif

// the fluid system
SET_PROP(ChannelTestTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
#if NONISOTHERMAL
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#else
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
#endif
};

// Set the grid type
SET_TYPE_PROP(ChannelTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ChannelTestTypeTag, Problem, Dumux::ChannelTestProblem<TypeTag> );

SET_BOOL_PROP(ChannelTestTypeTag, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(ChannelTestTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(ChannelTestTypeTag, EnableGridVolumeVariablesCache, true);


#if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(ChannelTestTypeTag, EnableInertiaTerms, true);
#else
SET_BOOL_PROP(ChannelTestTypeTag, EnableInertiaTerms, false);
#endif
}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 * \todo doc me!
 */
template <class TypeTag>
class ChannelTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    static constexpr auto dimWorld = GET_PROP_TYPE(TypeTag, GridView)::dimensionworld;

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

public:
    ChannelTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");

        using CellArray = std::array<unsigned int, dimWorld>;
        const CellArray numCells = getParam<CellArray>("Grid.Cells");
        cellSizeY_ = this->fvGridGeometry().bBoxMax()[1] / numCells[1];
    }

   /*!
     * \name Problem parameters
     */
    // \{


    bool shouldWriteRestartFile() const
    {
        return false;
    }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }
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

        if(isInlet(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
#if NONISOTHERMAL
            values.setDirichlet(Indices::temperatureIdx);
#endif
        }
        else if(isOutlet(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
#if NONISOTHERMAL
            values.setOutflow(Indices::energyBalanceIdx);
#endif
        }
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
#if NONISOTHERMAL
            values.setNeumann(Indices::energyBalanceIdx);
#endif
        }

        return values;
    }

   /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        if(isInlet(globalPos))
        {
#if NONISOTHERMAL
        // give the system some time so that the pressure can equilibrate, then start the injection of the hot liquid
        if(time() >= 200.0)
            values[Indices::temperatureIdx] = 293.15;
#endif
        }

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[Indices::pressureIdx] = 1.1e+5;
        values[Indices::velocityXIdx] = 0.0;
        values[Indices::velocityYIdx] = 0.0;

        if(isInlet(globalPos))
        {
            const auto& y = globalPos[1];
            const auto& yMax = this->fvGridGeometry().bBoxMax()[1];
            const auto& vMax = inletVelocity_;
            values[Indices::velocityXIdx] = 4 * vMax * y * (yMax - y)/(yMax * yMax);
        }

#if NONISOTHERMAL
        values[Indices::temperatureIdx] = 283.15;
#endif

        return values;
    }

    // \}
    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
        if(inletVelocity_ > eps_)
            timeLoop_->setCheckPoint({200.0, 210.0});
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

private:

    bool isInlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    bool isWall(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > eps_ || globalPos[0] < this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    bool isLowerLeftFace_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_ && globalPos[1] < (0.5*cellSizeY_ + eps_);
    }

    bool isUpperLeftFace_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_ && globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - 0.5 * cellSizeY_ - eps_;
    }

    Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
    Scalar cellSizeY_;
};
} //end namespace

#endif
