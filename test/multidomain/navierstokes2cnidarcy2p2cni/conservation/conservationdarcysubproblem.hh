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
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */
#ifndef DUMUX_CONSERVATION_DARCY_SUBPROBLEM_HH
#define DUMUX_CONSERVATION_DARCY_SUBPROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include "../properties.hh" // TODO
#include <dumux/porousmediumflow/implicit/problem.hh>
#include "conservationspatialparams.hh"

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

//#define ISOTHERMAL 0 // TODO NONISOTHERMAL ?!

namespace Dumux
{
template <class TypeTag>
class ConservationDarcySubProblem;

namespace Properties
{
NEW_TYPE_TAG(DarcySubProblem, INHERITS_FROM(CCTpfaModel, TwoPTwoCNI, ConservationSpatialParams));

// Set the problem property
SET_TYPE_PROP(DarcySubProblem, Problem, ConservationDarcySubProblem<TypeTag>);

// Set the grid type // TODO from darcytestproblem (1p)
#if ENABLE_3D
SET_TYPE_PROP(DarcySubProblem, Grid, Dune::YaspGrid<3>);
#else
SET_TYPE_PROP(DarcySubProblem, Grid, Dune::YaspGrid<2>);
#endif

// Disable gravity
SET_BOOL_PROP(DarcySubProblem, ProblemEnableGravity, false);

// choose pn and Sw as primary variables
SET_INT_PROP(DarcySubProblem, Formulation, TwoPTwoCFormulation::pnsw);

// the gas component balance (air) is replaced by the total mass balance
SET_INT_PROP(DarcySubProblem, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, Indices)::contiNEqIdx);

//// Set the grid type
//#if ENABLE_3D
//SET_TYPE_PROP(DarcySubProblem, Grid, Dune::YaspGrid< 3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
//#else
//SET_TYPE_PROP(DarcySubProblem, Grid, Dune::YaspGrid< 2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
//#endif

// Set the grid parameter group
SET_STRING_PROP(DarcySubProblem, GridParameterGroup, "DarcyGrid");

SET_BOOL_PROP(DarcySubProblem, EnableGlobalFVGeometryCache, true); // TODO default = false, but true needed for couplingmanager (access Darcy element)

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium. During
 *        buoyancy driven upward migration the gas passes a high
 *        temperature area.
 *
 */
template <class TypeTag>
class ConservationDarcySubProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
        temperatureIdx = Indices::temperatureIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    // property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ConservationDarcySubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        pressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, Pressure);
        switch_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, Switch);
        temperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, Temperature);
        initialPhasePresence_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, PorousMedium, InitialPhasePresence);
        name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);

        eps_ = 1e-6;

        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, DarcyGrid, UpperRight)[0];
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, DarcyGrid, UpperRight)[1];
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
    const std::string name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, NamePM); }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        ParentType::init();
    }

    // suppress output from DuMuX
    bool shouldWriteOutput() const
    {
        return true;
    }

    /*!
     * \brief Returns the source term
     *
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values = Scalar(0);

        return values;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllNeumann();

        double eps = 1.0e-5;
        if (globalPos[dim-1] < eps)
        {
          values.setDirichlet(pressureIdx);
          values.setDirichlet(switchIdx);
          values.setDirichlet(temperatureIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        initial_(values, globalPos);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values = Scalar(0);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        initial_(values, globalPos);

        return values;
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param globalPos The global position
     */
    int initialPhasePresenceAtPos(const GlobalPosition &globalPos) const
    {
        return initialPhasePresence_;
    }

    /*!
     * \brief Set the coupling manager
     * \param couplingManager The coupling manager
     *
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     *
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return globalPos[1] > bBoxMax_[1] - eps_; } // onUpperBoundary

private:
    /*!
     * \brief Internal method for the initial condition
     *        (reused for the dirichlet conditions!)
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values.setState(Indices::wPhaseOnly);
        values[pressureIdx] = pressure_;
        values[switchIdx] = switch_;
        values[temperatureIdx] = temperature_;
    }

    Scalar pressure_;
    Scalar switch_;
    Scalar temperature_;
    int initialPhasePresence_;
    std::string name_;

    Scalar eps_;

    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif // DUMUX_CONSERVATION_DARCY_SUBPROBLEM_HH
