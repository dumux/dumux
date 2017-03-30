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
 * \brief Channel flow test for the nonisothermal staggered grid (Navier-)Stokes model
 */
#ifndef DUMUX_CHANNEL_NI_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_NI_TEST_PROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggeredni/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
//#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
//#include <dumux/material/components/constant.hh>

namespace Dumux
{
template <class TypeTag>
class ChannelNITestProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<ChannelNITestProblem<TypeTag>>
    { static const bool value = false; };
}

namespace Properties
{
NEW_TYPE_TAG(ChannelNITestProblem, INHERITS_FROM(StaggeredModel, NavierStokesNI));

//SET_PROP(ChannelNITestProblem, Fluid)
//{
//private:
//    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//public:
//    typedef FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
//};

NEW_PROP_TAG(FluidSystem);

// Select the fluid system
SET_TYPE_PROP(ChannelNITestProblem, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)/*, SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)>, true*/>);

SET_PROP(ChannelNITestProblem, PhaseIdx)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    static constexpr int value = FluidSystem::wPhaseIdx;
};

// Set the grid type
SET_TYPE_PROP(ChannelNITestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ChannelNITestProblem, Problem, Dumux::ChannelNITestProblem<TypeTag> );

SET_BOOL_PROP(ChannelNITestProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(ChannelNITestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(ChannelNITestProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(ChannelNITestProblem, ProblemEnableGravity, true);

//#if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(ChannelNITestProblem, EnableInertiaTerms, true);
//#else
//SET_BOOL_PROP(ChannelNITestProblem, EnableInertiaTerms, false);
//#endif
}

/*!
 * \brief  Test problem for the one-phase model:
   \todo doc me!
 */
template <class TypeTag>
class ChannelNITestProblem : public NavierStokesProblem<TypeTag>
{
    typedef NavierStokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

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
        energyBalanceIdx = Indices::energyBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        temperatureIdx = Indices::temperatureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

public:
    ChannelNITestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        inletVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             Scalar,
                                             Problem,
                                             InletVelocity);
        FluidSystem::init();
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

    bool shouldWriteRestartFile() const
    {
        return false;
    }

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
        values.setDirichlet(energyBalanceIdx);

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
        BoundaryValues values;
        values[pressureIdx] = 1.1e+5;
        values[temperatureIdx]= 283.15; // TODO

        if(isInlet(globalPos))
        {
            values[velocityXIdx] = inletVelocity_;
            values[velocityYIdx] = 0.0;
            values[temperatureIdx] = 303.15; // TODO
        }
        else if(isWall(globalPos))
        {
            values[velocityXIdx] = 0.0;
            values[velocityYIdx] = 0.0;
        }
        else if(isOutlet(globalPos))
        {
            values[velocityXIdx] = 1.0;
            values[velocityYIdx]= 0.0;
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
        values[temperatureIdx] = 283.15;

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

    const Scalar eps_{1e-6};
    Scalar inletVelocity_;
    std::string name_;
};
} //end namespace

#endif
