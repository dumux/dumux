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


#include <dumux/linear/amgbackend.hh>


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
NEW_TYPE_TAG(ChannelTestProblem, INHERITS_FROM(StaggeredModel, NavierStokes));

SET_PROP(ChannelTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> > type;
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
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p -parameterFile test_box1p.input</tt> or
 * <tt>./test_cc1p -parameterFile test_cc1p.input</tt>
 *
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::YaspGrid<2> type;</tt> to
 * <tt>typedef Dune::YaspGrid<3> type;</tt> in the problem file
 * and use <tt>test_1p_3d.dgf</tt> in the parameter file.
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
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

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
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    CellCenterPrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return CellCenterPrimaryVariables(0);
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
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);

        // set a fixed pressure in one cell
        if (isOutlet(globalPos))
        {
            values.setDirichlet(massBalanceIdx);
            values.setOutflow(momentumBalanceIdx);
        }
        else
            values.setNeumann(massBalanceIdx);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    CellCenterPrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        CellCenterPrimaryVariables values(0);
        values[pressureIdx] = 1.1e+5;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    CellCenterPrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        return CellCenterPrimaryVariables(0);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    CellCenterPrimaryVariables initialCCValuesAtPos(const GlobalPosition &globalPos) const
    {
        CellCenterPrimaryVariables priVars(0);
        priVars[0] = 1.0e+5; //TODO: fix indices
        return priVars;
    }


    /*!
     * \brief Evaluate the initial value for a facet.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the initial values ought to be
     *            set (in global coordinates)
     * \param direction The direction index of the facets unit outer normal
     */
    GlobalPosition initialVelocityAtPos(const GlobalPosition &globalPos) const
    {
        GlobalPosition velocity;
        velocity[0] = 0.0;
        velocity[1] = 0.0;
        return velocity;

    }

     /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        facet.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the dirichlet condition ought to be
     *            set in global coordinates
     * \param direction The direction index of the facets unit outer normal
     */
    GlobalPosition dirichletVelocityAtPos(const GlobalPosition &pos) const
    {
        GlobalPosition velocity;

        if(isInlet(pos))
        {
            velocity[0] = inletVelocity_;
            velocity[1] = 0.0;
        }
        else if(isWall(pos))
        {
            velocity[0] = 0.0;
            velocity[1] = 0.0;
        }
        else if(isOutlet(pos))
        {
            velocity[0] = 1.0;
            velocity[1] = 0.0;
        }

        return velocity;
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
        return globalPos[1] < eps_ || globalPos[1] > this->bBoxMax()[1] - eps_;
    }


    Scalar eps_;
    Scalar inletVelocity_;
    std::string name_;
};
} //end namespace

#endif
