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
 * @brief  Definition of a simple Stokes problem
 */
#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggered/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

namespace Dumux
{

template <class TypeTag>
class StokesTestProblem;

//////////
// Specify the properties for the Stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesTestProblem, INHERITS_FROM(StaggeredModel, NavierStokes));

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
#else
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
#endif

// Set the problem property
SET_TYPE_PROP(StokesTestProblem, Problem, StokesTestProblem<TypeTag>);

SET_PROP(StokesTestProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> >;
};

SET_BOOL_PROP(StokesTestProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(StokesTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(StokesTestProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(StokesTestProblem, ProblemEnableGravity, false);

SET_BOOL_PROP(StokesTestProblem, EnableInertiaTerms, false);

// Set the grid parameter group
SET_STRING_PROP(StokesTestProblem, GridParameterGroup, "StokesGrid");

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);

}

/*!
 * \ingroup StaggeredStokesModel
 * \ingroup ImplicitTestProblems
 * \brief Stokes flow problem with water (H2O) flowing from the left to the right.
 *
 * The domain is sized 2m times 1m. The boundary conditions for the momentum balances
 * are set to Dirichlet with outflow on the right boundary. The mass balance has
 * outflow boundary conditions.
 *
 * This problem uses the \ref StaggeredStokesModel.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes</tt>
 */
template <class TypeTag>
class StokesTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using SubControlVolumeFace =  typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry) ;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

    using Fluid = typename GET_PROP_TYPE(TypeTag, Fluid);

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

public:
    StokesTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
{
        eps_ = 1e-6;
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        vIn_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, Velocity);
        pressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, FFPressure);

        bBoxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[0];
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[0];
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[1];
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[1];

        this->boundingBoxTree(); // TODO
}

    /*!
     * \name Problem parameters
     */
    // \{

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
        return 273.15 + 10.0; // 10°C
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // inflow left, outflow right
        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);

        // set a fixed pressure in one cell
        if (onRightBoundary_(globalPos))
        {
            values.setDirichlet(massBalanceIdx);
            values.setOutflow(momentumBalanceIdx);
        }
        else
            values.setOutflow(massBalanceIdx);

        if (onCouplingInterface(globalPos))
        {
            values.setDirichlet(massBalanceIdx);
            values.setCouplingNeumann(momentumBalanceIdx);
        }

        return values;
    }

    //! \copydoc ImplicitProblem::dirichlet()
    BoundaryValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.dofPosition();

        BoundaryValues values(0.0);

        values = initialAtPos(globalPos);

        return values;
    }

    //! \copydoc ImplicitProblem::neumann()
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

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);
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

        values[pressureIdx] = pressure_;
        // left inflow
        values[velocityXIdx] = vIn_*(globalPos[1] - bBoxMin_[1])*(bBoxMax_[1] - globalPos[1])
                                 / (0.25*(bBoxMax_[1] - bBoxMin_[1])*(bBoxMax_[1] - bBoxMin_[1]));
        return values;
    }

    // \}

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
    {return onLowerBoundary_(globalPos); }

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
    Scalar vIn_;
    Scalar pressure_;
    std::string name_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif
