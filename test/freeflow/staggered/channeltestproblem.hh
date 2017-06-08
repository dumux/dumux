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
 *
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model
 */
#ifndef DUMUX_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_TEST_PROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggered/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>

namespace Dumux
{
template <class TypeTag>
class ChannelTestProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<ChannelTestProblem<TypeTag>>
    { static const bool value = false; };
}

namespace Properties
{
#if !NONISOTHERMAL
NEW_TYPE_TAG(ChannelTestProblem, INHERITS_FROM(StaggeredModel, NavierStokes));
#else
NEW_TYPE_TAG(ChannelTestProblem, INHERITS_FROM(StaggeredModel, NavierStokesNI));
#endif

SET_PROP(ChannelTestProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
#if NONISOTHERMAL
    using type = FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > ;
#else
    using type = FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> > ;
#endif
};

// Set the grid type
SET_TYPE_PROP(ChannelTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ChannelTestProblem, Problem, Dumux::ChannelTestProblem<TypeTag> );

SET_BOOL_PROP(ChannelTestProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(ChannelTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(ChannelTestProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(ChannelTestProblem, ProblemEnableGravity, true);

#if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(ChannelTestProblem, EnableInertiaTerms, true);
#else
SET_BOOL_PROP(ChannelTestProblem, EnableInertiaTerms, false);
#endif
}

/*!
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel:
   \todo doc me!
 */
template <class TypeTag>
class ChannelTestProblem : public NavierStokesProblem<TypeTag>
{
    typedef NavierStokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyBalanceIdx = Indices::energyBalanceIdx,
#endif
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

public:
    ChannelTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView), eps_(1e-6)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        inletVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             Scalar,
                                             Problem,
                                             InletVelocity);

#if NONISOTHERMAL
    if(inletVelocity_ > eps_)
        this->timeManager().startNextEpisode(200.0);
#endif
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    {
        return name_;
    }

#if NONISOTHERMAL
    void episodeEnd()
    {
        if(inletVelocity_ > eps_)
        {
            this->timeManager().startNextEpisode(50.0);
            this->timeManager().setTimeStepSize(10.0);
        }
    }
#endif

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
    SourceValues sourceAtPos(const GlobalPosition &globalPos) const
    {
        return SourceValues(0.0);
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

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);

#if NONISOTHERMAL
        if(isInlet(globalPos))
            values.setDirichlet(energyBalanceIdx);
        else
            values.setOutflow(energyBalanceIdx);
#endif

        // set a fixed pressure in one cell
        if (isOutlet(globalPos))
        {
            values.setDirichlet(massBalanceIdx);
            values.setOutflow(momentumBalanceIdx);
        }
        else
            values.setOutflow(massBalanceIdx);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    BoundaryValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryValues values = initialAtPos(globalPos);

        if(isInlet(globalPos))
        {
            values[velocityXIdx] = inletVelocity_;
#if NONISOTHERMAL
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();

        // give the system some time so that the pressure can equilibrate, then start the injection of the hot liquid
        if(time > 200.0)
            values[temperatureIdx] = 293.15;
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
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values;
        values[pressureIdx] = 1.1e+5;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

#if NONISOTHERMAL
        values[temperatureIdx] = 283.15;
#endif

        return values;
    }

    // \}

private:

    bool isInlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool isWall(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > eps_ || globalPos[0] < this->bBoxMax()[0] - eps_;
    }

    Scalar eps_;
    Scalar inletVelocity_;
    std::string name_;
};
} //end namespace

#endif
