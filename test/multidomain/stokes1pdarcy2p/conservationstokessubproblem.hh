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

#ifndef DUMUX_CONSERVATION_STOKES_SUBPROBLEM_HH
#define DUMUX_CONSERVATION_STOKES_SUBPROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggered/model.hh>
#include <dumux/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
// test
#include <dumux/material/components/air.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/components/constant.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

namespace Dumux
{
template <class TypeTag>
class ConservationStokesSubProblem;

//////////
// Specify the properties for the Stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesSubProblem, INHERITS_FROM(StaggeredModel, NavierStokes));

// Set the problem property
SET_TYPE_PROP(StokesSubProblem, Problem, ConservationStokesSubProblem<TypeTag>);

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(StokesSubProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
#else
SET_TYPE_PROP(StokesSubProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
#endif

SET_TYPE_PROP(StokesSubProblem, FluidSystem, FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// set phase index (air)
SET_INT_PROP(StokesSubProblem, PhaseIdx, GET_PROP_TYPE(TypeTag, FluidSystem)::nPhaseIdx);

SET_BOOL_PROP(StokesSubProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(StokesSubProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(StokesSubProblem, EnableGlobalVolumeVariablesCache, true);

// Use complete Navier-Stokes equation
SET_BOOL_PROP(StokesSubProblem, EnableInertiaTerms, false);

// Disable gravity field
SET_BOOL_PROP(StokesSubProblem, ProblemEnableGravity, false);

// TODO define in global problem?
SET_BOOL_PROP(StokesSubProblem, UseMoles, false);

// Set the grid parameter group
SET_STRING_PROP(StokesSubProblem, GridParameterGroup, "StokesGrid");

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);
}

template <class TypeTag>
class ConservationStokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum { dim = GridView::dimension,
        dimWorld = GridView::dimensionworld};

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry) ;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);


public:
    ConservationStokesSubProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        velocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, Velocity);
        pressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, Pressure);
        temperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, Temperature);
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);

        bBoxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[0];
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[0];
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[1];
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[1];

        this->boundingBoxTree(); // TODO
    }

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     *
     */
    const std::string name() const
    {
        return name_+"_stokes";
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The global position at which the temperature is set
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return temperature_;
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        if (onCouplingInterface(globalPos))
        {
            values.setCouplingNeumann(momentumBalanceIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element
     * \param scvf The sub control volume face
     */
    BoundaryValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.dofPosition();

        BoundaryValues values(0.0);
        values = initialAtPos(globalPos);

        return values;
    }

    BoundaryValues neumann(const Element& element,
            const FVElementGeometry& fvGeometry,
            const ElementVolumeVariables& elemVolVars,
            const SubControlVolumeFace& scvf) const
    {
        BoundaryValues values(0.0);

        if(onCouplingInterface(scvf.center()))
        {
            values[velocityYIdx] = couplingManager().darcyData().momentumCouplingCondition(scvf);
        }
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values(0.0);

        // pressure
        values[pressureIdx] = pressure_;

        // velocity
        if (dim == 1)
            values[velocityXIdx] = velocity_;
        else // dim == 2
            values[velocityYIdx] = 4.0 * velocity_ * globalPos[0] * (bBoxMax_[0] - globalPos[0])
                                   / (bBoxMax_[0] - bBoxMin_[0])
                                   / (bBoxMax_[0] - bBoxMin_[0]);

        return values;
    }

    /*!
     * \brief Set the coupling manager
     *
     * \param couplingManager The coupling manager
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    /*!
     * \brief Check if on coupling interface
     *
     * \param globalPos The global position
     *
     * Returns true if globalPos is on coupling interface
     * (here: lower boundary of Stokes domain)
     */
    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {
        if (dim == 1)
            return globalPos[0] < bBoxMin_[0] + eps_;
        else // dim == 2
            return globalPos[dim-1] < bBoxMin_[dim-1] + eps_;
    }

private:

    std::string name_;
    static constexpr Scalar eps_ = 1e-6;

    Scalar velocity_;
    Scalar pressure_;
    Scalar temperature_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif // DUMUX_CONSERVATION_STOKES_SUBPROBLEM_HH
