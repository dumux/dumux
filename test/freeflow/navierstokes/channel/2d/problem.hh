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
// Create new type tags
namespace TTag {
#if !NONISOTHERMAL
struct ChannelTest { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
#else
struct ChannelTest { using InheritsFrom = std::tuple<NavierStokesNI, StaggeredFreeFlowModel>; };
#endif
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
#if NONISOTHERMAL
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
#else
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
#endif
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelTest> { using type = Dumux::ChannelTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::ChannelTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelTest> { static constexpr bool value = true; };
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

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto dimWorld = GetPropType<TypeTag, Properties::GridView>::dimensionworld;

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

public:
    ChannelTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
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
            values.setOutflow(Indices::energyEqIdx);
#endif
        }
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
#if NONISOTHERMAL
            values.setNeumann(Indices::energyEqIdx);
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
            values[Indices::velocityXIdx] = inletVelocity_;
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

    Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
};
} //end namespace

#endif
