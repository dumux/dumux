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
/**
 * @file
 * @brief  Definition of a simple Darcy problem
 */
#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include "1pspatialparams.hh"
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

namespace Dumux
{
template <class TypeTag>
class DarcyTestProblem;

namespace Properties
{
NEW_TYPE_TAG(DarcyTestProblem, INHERITS_FROM(CCTpfaModel, OneP, OnePSpatialParams));

// Set the problem property
SET_TYPE_PROP(DarcyTestProblem, Problem, Dumux::DarcyTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(DarcyTestProblem, SpatialParams, OnePSpatialParams<TypeTag>);

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(DarcyTestProblem, Grid, Dune::YaspGrid<3>);
#else
SET_TYPE_PROP(DarcyTestProblem, Grid, Dune::YaspGrid<2>);
#endif

SET_PROP(DarcyTestProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> >;
};

// Enable gravity
SET_BOOL_PROP(DarcyTestProblem, ProblemEnableGravity, false);

// Set the grid parameter group
SET_STRING_PROP(DarcyTestProblem, GridParameterGroup, "DarcyGrid");

SET_BOOL_PROP(DarcyTestProblem, EnableGlobalFVGeometryCache, true); // TODO default = false, but true needed for couplingmanager (access Darcy element)

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
}

template <class TypeTag>
class DarcyTestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables) ;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry) ;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        // primary variable indices
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

    enum { dofCodim = 0 };

public:
    DarcyTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), eps_(1e-7)
{
        // get some parameters from the input file
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        bBoxMin_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, DarcyGrid, LowerLeft);
        bBoxMax_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, DarcyGrid, UpperRight);
}

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    {
        return name_+"_darcy";
    }

    bool shouldWriteOutput() const //define output
    {
        // write initial conditions
        return (this->timeManager().time() < 0.0);
    }

    void preTimeStep()
    {
    }

    /*!
     * \brief Called at the end of each time step
     */
    void postTimeStep()
    {
    }


    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if(onCouplingInterface(globalPos))
            values.setAllCouplingNeumann();

        return values;
    }


    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex (pore body) for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 1.0e5;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex (pore body) for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolvars,
                             const SubControlVolumeFace& scvf) const
    {
        return PrimaryVariables(0.0);
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        return PrimaryVariables(0.0);
    }
    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Element &element) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 1.0e5;
        return values;
    }


    // \}

    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    //! Get the coupling manager
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    /*!
     * \brief Return the coupling boundary
     *
     * \param globalPos The global position
     */
    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return onUpperBoundary_(globalPos); }

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    Scalar eps_;
    std::string name_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif
